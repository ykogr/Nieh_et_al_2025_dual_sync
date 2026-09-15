#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Video annotator (grouped version + undo last entry + real-time speed + seek backward/forward)
- Playback speed is calibrated to the video FPS, so 1 second is a real second.
- Supports "seek forward/backward 1 s": A (backward), D (forward)
- A "reach" starts a new row; subsequent "chews" are appended to that row.
- Supports "delete the last entry" (Undo): Backspace or 'i'
- Output: one CSV per video; columns: trial number, chew 1, chew 2, ...

Controls:
 [F] reach (add a new row)
 [J] chew (append to the current row)
 [I] or [Backspace] delete the last entry (can be pressed repeatedly)
 [A] back 1 s
 [D] forward 1 s
 [Space] pause/resume playback
 [N] next video
 [Q] quit the program (output is written automatically)
"""

import os
import sys
import glob
import csv
import argparse
from typing import List, Dict, Any
from datetime import timedelta
import cv2

VIDEO_EXTS = ['*.mp4', '*.mov', '*.avi', '*.mkv', '*.m4v', '*.wmv']

ACTION_REACH = ord('f')  # reach
ACTION_CHEW = ord('j')  # chew
PAUSE_KEY = 32  # Space
NEXT_KEY = ord('n')
QUIT_KEY = ord('q')
UNDO_KEYS = {8, 127, ord('i')}  # Backspace (8/127 depending on the system) or 'i'
BACKWARD_KEY = ord('a')  # back 1 s
FORWARD_KEY = ord('d')  # forward 1 s


def find_videos(folder: str) -> List[str]:
    paths = []
    for ext in VIDEO_EXTS:
        paths.extend(glob.glob(os.path.join(folder, ext)))
    paths.sort()
    return paths


def ms_to_ts(ms: float) -> str:
    seconds = ms / 1000.0
    td = timedelta(seconds=seconds)
    total_seconds = int(td.total_seconds())
    h = total_seconds // 3600
    m = (total_seconds % 3600) // 60
    s = total_seconds % 60
    millis = int(round((seconds - int(seconds)) * 1000))
    return f"{h:02d}:{m:02d}:{s:02d}.{millis:03d}"


def draw_hud(frame, filename: str, current_ms: float, reach_count: int, last_chew_count: int):
    lines = [
        f"Filename: {os.path.basename(filename)}",
        f"Time: {ms_to_ts(current_ms)}",
        f"Now_row(get chip): {reach_count} , record_chew: {last_chew_count} ",
        "Op: [f]intake [j]chew [i]undo [a]/[d] -/+1s [Space]stop [N]next [Q]end",
    ]
    y = 30
    for line in lines:
        cv2.putText(frame, line, (15, y), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 255), 2, cv2.LINE_AA)
        y += 30
    return frame

def play_and_collect(video_path: str):
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        print(f"Cannot open video: {video_path}")
        return [], False

    # Get the video FPS to compute the real playback delay (fixes playback being too fast)
    fps = cap.get(cv2.CAP_PROP_FPS)
    if not fps or fps <= 0:
        fps = 30.0  # safe default
    frame_delay = int(1000 / fps)

    # Get the total video length (ms) so seeking cannot go out of range
    total_frames = cap.get(cv2.CAP_PROP_FRAME_COUNT)
    total_ms = (total_frames / fps) * 1000 if fps else float('inf')

    cv2.namedWindow('Marker', cv2.WINDOW_NORMAL)
    cv2.resizeWindow('Marker', 1280, 720)

    rows: List[List[str]] = []  # dynamic rows
    history: List[Dict[str, Any]] = []  # action log, used for Undo
    paused = False

    # Pre-read the first frame
    ret, frame = cap.read()

    def flash_text(text: str, img_frame):
        if img_frame is not None:
            img = img_frame.copy()
            cv2.putText(img, text, (15, 180), cv2.FONT_HERSHEY_SIMPLEX, 1.0, (0, 200, 255), 3, cv2.LINE_AA)
            cv2.imshow('Marker', img)
            cv2.waitKey(220)

    while True:
        current_ms = cap.get(cv2.CAP_PROP_POS_MSEC)

        reach_count = len(rows)
        last_chew_count = len(rows[-1]) - 1 if rows else 0

        if frame is not None:
            hud = draw_hud(frame.copy(), video_path, current_ms, reach_count, last_chew_count)
            cv2.imshow('Marker', hud)

        # Wait time is set from the FPS; while paused, poll for keys at a 50 ms refresh rate
        key = cv2.waitKey(frame_delay if not paused else 50) & 0xFF

        if key != 255:  # a key was pressed
            if key == ACTION_REACH:
                rows.append([str(len(rows) + 1)])
                history.append({'type': 'reach', 'row': len(rows) - 1})
                flash_text("Marked: new_chip", frame)

            elif key == ACTION_CHEW:
                if not rows:
                    flash_text("Press [f] to build New_chip first", frame)
                else:
                    ts = f"{current_ms / 1000:.3f}"  # record to three decimal places (ms)
                    rows[-1].append(ts)
                    history.append({'type': 'chew', 'row': len(rows) - 1, 'ts': ts})
                    flash_text("Marked: chew", frame)

            elif key in UNDO_KEYS:
                if not history:
                    flash_text("No items can be deleted", frame)
                else:
                    last = history.pop()
                    if last['type'] == 'chew':
                        r = last['row']
                        if len(rows[r]) > 1 and rows[r][-1] == last.get('ts'):
                            rows[r].pop()
                        elif len(rows[r]) > 1:
                            rows[r].pop()
                        flash_text("Deleted: last chew", frame)
                    elif last['type'] == 'reach':
                        r = last['row']
                        if r == len(rows) - 1 and len(rows[r]) == 1:
                            rows.pop()
                            flash_text("Deleted: last new_chip", frame)
                        else:
                            flash_text("Cannot delete (has data in row)", frame)

            elif key == BACKWARD_KEY:  # A key: back 1 s
                new_ms = max(0, current_ms - 1000)
                cap.set(cv2.CAP_PROP_POS_MSEC, new_ms)
                ret, frame = cap.read()  # refresh the frame immediately
                flash_text("-1 Sec", frame)

            elif key == FORWARD_KEY:  # D key: forward 1 s
                new_ms = min(total_ms, current_ms + 1000)
                cap.set(cv2.CAP_PROP_POS_MSEC, new_ms)
                ret, frame = cap.read()  # refresh the frame immediately
                flash_text("+1 Sec", frame)

            elif key == PAUSE_KEY:  # * pause function restored
                paused = not paused

            elif key == NEXT_KEY:
                # Press N to go to the next video (without ending the program)
                cap.release()
                return rows, False

            elif key == QUIT_KEY:
                # Press Q to quit: first return the annotated data, then tell main() to exit
                cap.release()
                return rows, True

        # If not paused, keep reading the next frame
        if not paused and key == 255:
            ret, frame = cap.read()
            if not ret:
                break  # video finished playing

    cap.release()
    # * Critical fix: when the video ends naturally, rows must be returned together with False ("do not quit")
    return rows, False

def export_rows(video_path: str, rows: List[List[str]], out_dir: str):
    if not rows:
        return None
    max_chews = max((len(r) - 1) for r in rows)
    headers = ['Trial'] + [f'Chew_{i}' for i in range(1, max_chews + 1)]
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(video_path))[0]
    csv_path = os.path.join(out_dir, f'{base}_B.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for r in rows:
            padded = r + [''] * (1 + max_chews - len(r))
            writer.writerow(padded)
    return csv_path


def main():
    ap = argparse.ArgumentParser(description='Video annotator (grouped version + Undo + speed calibration)')
    ap.add_argument('folder', nargs='?', default='.', help='video folder (default: current folder)')
    ap.add_argument('--outdir', default='csv_outputs', help='output folder (one CSV per video)')
    args = ap.parse_args()

    videos = find_videos(args.folder)
    if not videos:
        print(f"[ERROR] No videos found in '{args.folder}' (supported: {', '.join(VIDEO_EXTS)})")
        sys.exit(1)

    print("The following videos will be played in order:")
    for v in videos:
        print(" -", os.path.basename(v))

    exported = []
    try:
        for v in videos:
            # Receive the annotated data (rows) and the quit signal (quit_flag)
            rows, quit_flag = play_and_collect(v)

            # Whenever there is annotated data, always save it first!
            path = export_rows(v, rows, args.outdir)
            if path:
                exported.append(path)
                print(f"[DONE] Written: {path}")

            # After saving, if the quit signal (Q) was received, stop processing the remaining videos
            if quit_flag:
                print("\nQuit command (Q) received; stopping processing of the remaining videos...")
                break

    except KeyboardInterrupt:
        # Safety net for a forced Ctrl+C in the terminal
        print("\nForced interrupt received; writing the files collected so far...")
    finally:
        cv2.destroyAllWindows()

    if exported:
        print("\nOutput list:")
        for p in exported:
            print(" -", p)
    else:
        print("No data to output.")

if __name__ == '__main__':
    main()