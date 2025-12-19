#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Video annotator
- Press "reach" to start a new row; subsequent "chew" marks are appended to that row.
- Supports undo of the most recent action: Backspace or 'i'
- Output: one CSV per video; columns: trial index, chew1, chew2, ...

Controls:
 [f] Reach (start a new row)
 [j] Chew (append to the current row)
 [Backspace]/[i] Undo the most recent action (can be repeated)
 [Space] Pause/Resume
 [N] Next video
 [Q] Quit (auto-export)
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

ACTION_REACH = ord('f')   # Reach
ACTION_CHEW  = ord('j')   # Chew
PAUSE_KEY    = 32         # Space
NEXT_KEY     = ord('n')
QUIT_KEY     = ord('q')
UNDO_KEYS    = {8, 127, ord('i')}  # Backspace (8/127 depending on OS) or 'i'

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
        "Operation: [f]intake  [j]chew  [i]undo  [Space]stop  [N]next video  [Q]end",
    ]
    y = 30
    for line in lines:
        cv2.putText(frame, line, (15, y), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255,255,255), 2, cv2.LINE_AA)
        y += 30
    return frame

def play_and_collect(video_path: str) -> List[List[str]]:
    """
    Returns rows: each row is [trial_index, chew1, chew2, ...] (timestamps as strings).
    Also supports Undo: reverts the most recent addition (reach or chew).
    """
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        print(f"Cannot open video: {video_path}")
        return []

    cv2.namedWindow('Marker', cv2.WINDOW_NORMAL)
    cv2.resizeWindow('Marker', 1280, 720)

    rows: List[List[str]] = []   # Dynamic rows
    history: List[Dict[str, Any]] = []  # Action history for Undo
    paused = False

    def flash_text(text: str, frame):
        if frame is not None:
            img = frame.copy()
            cv2.putText(img, text, (15, 150), cv2.FONT_HERSHEY_SIMPLEX, 1.0, (0,200,255), 3, cv2.LINE_AA)
            cv2.imshow('Marker', img)
            cv2.waitKey(220)

    while True:
        if not paused:
            ret, frame = cap.read()
            if not ret:
                break
        current_ms = cap.get(cv2.CAP_PROP_POS_MSEC)

        reach_count = len(rows)
        last_chew_count = len(rows[-1]) - 1 if rows else 0

        if frame is not None:
            hud = draw_hud(frame.copy(), video_path, current_ms, reach_count, last_chew_count)
            cv2.imshow('Marker', hud)

        key = cv2.waitKey(1 if not paused else 50) & 0xFF
        if key == 255:
            continue

        if key == ACTION_REACH:
            rows.append([str(len(rows)+1)])
            history.append({'type': 'reach', 'row': len(rows)-1})
            flash_text("Marked: new chip）", frame)

        elif key == ACTION_CHEW:
            if not rows:
                flash_text("Please mark a new chip first", frame)
                continue
            ts = f"{current_ms/1000:.1f}"
            rows[-1].append(ts)
            history.append({'type': 'chew', 'row': len(rows)-1, 'ts': ts})
            flash_text("Marked: chew", frame)

        elif key in UNDO_KEYS:
            if not history:
                flash_text("Nothing to undo", frame)
                continue
            last = history.pop()
            if last['type'] == 'chew':
                r = last['row']
                # Safety check: delete only the last timestamp and confirm it matches
                if len(rows[r]) > 1 and rows[r][-1] == last.get('ts'):
                    rows[r].pop()
                elif len(rows[r]) > 1:
                    rows[r].pop()  # If it doesn't match, still pop the last one
                flash_text("Deleted: last chew", frame)
            elif last['type'] == 'reach':
                r = last['row']
                # Only delete the whole row if it still contains only the trial index
                if r == len(rows)-1 and len(rows[r]) == 1:
                    rows.pop()
                    flash_text("Undone: last new chip", frame)
                else:
                    # Should not happen: after a reach mark you would normally add chews, which would then appear in history
                    flash_text("Cannot undo reach (row already has data)", frame)

        elif key == PAUSE_KEY:
            paused = not paused

        elif key == NEXT_KEY:
            break

        elif key == QUIT_KEY:
            cap.release()
            raise KeyboardInterrupt

    cap.release()
    return rows

def export_rows(video_path: str, rows: List[List[str]], out_dir: str):
    if not rows:
        return None
    max_chews = max((len(r) - 1) for r in rows)
    headers = ["Trial"] + [f"Chew{i}" for i in range(1, max_chews + 1)]
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(video_path))[0]
    csv_path = os.path.join(out_dir, f"{base}_B.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for r in rows:
            padded = r + [""] * (1 + max_chews - len(r))
            writer.writerow(padded)
    return csv_path


def main():
    ap = argparse.ArgumentParser(description="Video annotator (grouped + Undo)")
    ap.add_argument("folder", nargs="?", default=".", help="Video folder (default: current directory)")
    ap.add_argument("--outdir", default="csv_outputs", help="Output folder (one CSV per video)")
    args = ap.parse_args()

    videos = find_videos(args.folder)
    if not videos:
        print(f"[ERROR] No videos found in '{args.folder}' (supported: {', '.join(VIDEO_EXTS)})")
        sys.exit(1)

    print("Videos will be played in the following order:")
    for v in videos:
        print(" -", os.path.basename(v))

    exported = []
    try:
        for v in videos:
            rows = play_and_collect(v)
            path = export_rows(v, rows, args.outdir)
            if path:
                exported.append(path)
                print(f"[DONE] Exported: {path}")
    except KeyboardInterrupt:
        print("\nReceived quit signal. Preparing to finalize exported files...")
    finally:
        cv2.destroyAllWindows()

    if exported:
        print("\nExport list:")
        for p in exported:
            print(" -", p)
    else:
        print("No data to export.")


if __name__ == "__main__":
    main()

