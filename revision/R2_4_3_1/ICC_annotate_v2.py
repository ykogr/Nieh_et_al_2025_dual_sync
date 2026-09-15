#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ICC_annotate_v2.py — annotation tool for the second coder (GUI-button version)

Additions over the original ICC_annotate.py:
- Clickable control buttons (no keyboard needed; keys still work as well)
- Slow playback: x0.25 / x0.5 / x1
- Frame step forward/backward while paused (< / >)
- Fine ±0.2 s nudge (±1 s as before)

Controls (button or key):
  INTAKE [F] reach (start a new row)   CHEW [J] chew (one jaw closure = 1 event)
  UNDO [I/Backspace] undo the last entry   PAUSE [Space] pause/resume
  -1s [A] / +1s [D]   -0.2s [Z] / +0.2s [C]   <FRAME [,] / FRAME> [.]
  x0.25 [1] / x0.5 [2] / x1 [3]   NEXT [N] next video   QUIT [Q] save and quit

Output: one <video name>_B.csv per video (columns: Trial, Chew_1, Chew_2, ...)
Usage: python ICC_annotate_v2.py <video folder> --outdir <output folder>
"""

import os
import sys
import glob
import csv
import time
import argparse
from typing import List, Dict, Any
from datetime import timedelta
import cv2
import numpy as np

VIDEO_EXTS = ['*.mp4', '*.mov', '*.avi', '*.mkv', '*.m4v', '*.wmv']

ACTION_REACH = ord('f')
ACTION_CHEW = ord('j')
PAUSE_KEY = 32  # Space
NEXT_KEY = ord('n')
QUIT_KEY = ord('q')
UNDO_KEY = ord('i')
UNDO_KEYS = {8, 127, UNDO_KEY}
BACKWARD_KEY = ord('a')
FORWARD_KEY = ord('d')
SMALL_BACK_KEY = ord('z')
SMALL_FWD_KEY = ord('c')
FRAME_BACK_KEY = ord(',')
FRAME_FWD_KEY = ord('.')
SPEED_KEYS = {ord('1'): 0.25, ord('2'): 0.5, ord('3'): 1.0}

DISP_W = 1280          # fixed width of the display canvas (keeps button coordinates stable)
PANEL_H = 150          # height of the button panel
WINDOW = 'Marker'


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


def build_buttons(panel_top: int):
    """Button definitions: (x1, y1, x2, y2, label, keycode). y coordinates are relative to the whole canvas."""
    btns = []
    # --- top row (large annotation buttons) ---
    row1_y1, row1_y2 = panel_top + 6, panel_top + 70
    row1 = [("INTAKE [F]", ACTION_REACH, 240),
            ("CHEW [J]", ACTION_CHEW, 420),
            ("UNDO [I]", UNDO_KEY, 240),
            ("PAUSE [Sp]", PAUSE_KEY, 240)]
    x = 10
    for label, code, w in row1:
        btns.append((x, row1_y1, x + w, row1_y2, label, code))
        x += w + 10
    # --- bottom row (seek, speed, control) ---
    row2_y1, row2_y2 = panel_top + 78, panel_top + 142
    row2 = [("-1s", BACKWARD_KEY), ("-0.2", SMALL_BACK_KEY),
            ("<FR", FRAME_BACK_KEY), ("FR>", FRAME_FWD_KEY),
            ("+0.2", SMALL_FWD_KEY), ("+1s", FORWARD_KEY),
            ("x.25", ord('1')), ("x.5", ord('2')), ("x1", ord('3')),
            ("NEXT", NEXT_KEY), ("QUIT", QUIT_KEY)]
    n = len(row2)
    gap = 8
    w = (DISP_W - 20 - gap * (n - 1)) // n
    x = 10
    for label, code in row2:
        btns.append((x, row2_y1, x + w, row2_y2, label, code))
        x += w + gap
    return btns


def draw_panel(canvas, buttons, speed: float, paused: bool):
    for (x1, y1, x2, y2, label, code) in buttons:
        hot = (code == ACTION_CHEW)
        sel = (code in SPEED_KEYS and SPEED_KEYS[code] == speed) or \
              (code == PAUSE_KEY and paused)
        fill = (60, 130, 60) if hot else ((90, 90, 40) if sel else (70, 70, 70))
        cv2.rectangle(canvas, (x1, y1), (x2, y2), fill, -1)
        cv2.rectangle(canvas, (x1, y1), (x2, y2), (200, 200, 200), 1)
        fs = 0.85 if hot else 0.55
        th = 2 if hot else 1
        (tw, tht), _ = cv2.getTextSize(label, cv2.FONT_HERSHEY_SIMPLEX, fs, th)
        tx = x1 + max(4, ((x2 - x1) - tw) // 2)
        ty = y1 + ((y2 - y1) + tht) // 2
        cv2.putText(canvas, label, (tx, ty), cv2.FONT_HERSHEY_SIMPLEX, fs,
                    (255, 255, 255), th, cv2.LINE_AA)


def draw_hud(canvas, filename, current_ms, reach_count, last_chew_count,
             speed, paused, flash_msg, last_chew_ts=None):
    state = f"Speed x{speed:g}" + ("  [PAUSED]" if paused else "")
    last_str = f"   LAST CHEW: {float(last_chew_ts):.1f}s" if last_chew_ts else "   LAST CHEW: ---"
    lines = [
        f"{os.path.basename(filename)}   Time: {ms_to_ts(current_ms)}   {state}",
        f"Now_row(get chip): {reach_count} , record_chew: {last_chew_count}{last_str}",
    ]
    y = 30
    for line in lines:
        cv2.putText(canvas, line, (15, y), cv2.FONT_HERSHEY_SIMPLEX, 0.7,
                    (255, 255, 255), 2, cv2.LINE_AA)
        y += 30
    if flash_msg:
        cv2.putText(canvas, flash_msg, (15, y + 20), cv2.FONT_HERSHEY_SIMPLEX,
                    1.0, (0, 200, 255), 3, cv2.LINE_AA)


def load_existing(video_path: str, out_dir: str) -> List[List[str]]:
    """If an autosaved CSV with the same name exists, load it so annotation can resume from there."""
    base = os.path.splitext(os.path.basename(video_path))[0]
    csv_path = os.path.join(out_dir, f'{base}_B.csv')
    if not os.path.exists(csv_path):
        return []
    rows: List[List[str]] = []
    try:
        with open(csv_path, newline='', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            for r in reader:
                cells = [c for c in r if c.strip() != '']
                if cells:
                    rows.append(cells)
    except Exception:
        return []
    return rows


def play_and_collect(video_path: str, autosave_dir: str = None):
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        print(f"Cannot open video: {video_path}")
        return [], False

    fps = cap.get(cv2.CAP_PROP_FPS)
    if not fps or fps <= 0:
        fps = 30.0
    total_frames = cap.get(cv2.CAP_PROP_FRAME_COUNT)
    total_ms = (total_frames / fps) * 1000 if fps else float('inf')

    ret, frame = cap.read()
    if frame is None:
        print(f"Cannot read first frame: {video_path}")
        cap.release()
        return [], False

    h0, w0 = frame.shape[:2]
    disp_h = int(h0 * DISP_W / w0)
    buttons = build_buttons(panel_top=disp_h)
    canvas_h = disp_h + PANEL_H

    cv2.namedWindow(WINDOW, cv2.WINDOW_AUTOSIZE)  # fixed size = stable click coordinates

    clicks: List[int] = []

    def on_mouse(event, mx, my, flags, param):
        if event == cv2.EVENT_LBUTTONDOWN:
            for (x1, y1, x2, y2, label, code) in buttons:
                if x1 <= mx <= x2 and y1 <= my <= y2:
                    clicks.append(code)
                    return

    cv2.setMouseCallback(WINDOW, on_mouse)

    rows: List[List[str]] = load_existing(video_path, autosave_dir) if autosave_dir else []
    history: List[Dict[str, Any]] = []
    paused = False
    speed = 1.0
    flash_msg, flash_until = "", 0.0
    if rows:
        print(f"Resumed from existing CSV: {len(rows)} rows loaded")

    def flash(text):
        nonlocal flash_msg, flash_until
        flash_msg, flash_until = text, time.time() + 0.6

    while True:
        current_ms = cap.get(cv2.CAP_PROP_POS_MSEC)
        reach_count = len(rows)
        last_chew_count = len(rows[-1]) - 1 if rows else 0

        disp = cv2.resize(frame, (DISP_W, disp_h))
        canvas = np.zeros((canvas_h, DISP_W, 3), dtype=np.uint8)
        canvas[:disp_h] = disp
        if time.time() > flash_until:
            flash_msg = ""
        last_chew_ts = None
        for r in reversed(rows):
            if len(r) > 1:
                last_chew_ts = r[-1]
                break
        draw_hud(canvas, video_path, current_ms, reach_count, last_chew_count,
                 speed, paused, flash_msg, last_chew_ts)
        draw_panel(canvas, buttons, speed, paused)
        cv2.imshow(WINDOW, canvas)

        frame_delay = max(1, int(1000 / fps / speed))
        key = cv2.waitKey(frame_delay if not paused else 50) & 0xFF
        if key == 255 and clicks:
            key = clicks.pop(0)
        marked = key in (ACTION_REACH, ACTION_CHEW) or key in UNDO_KEYS

        if key != 255:
            if key == ACTION_REACH:
                rows.append([str(len(rows) + 1)])
                history.append({'type': 'reach', 'row': len(rows) - 1})
                flash("Marked: new_chip")

            elif key == ACTION_CHEW:
                if not rows:
                    flash("Press INTAKE first")
                else:
                    ts = f"{current_ms / 1000:.3f}"
                    rows[-1].append(ts)
                    history.append({'type': 'chew', 'row': len(rows) - 1, 'ts': ts})
                    flash(f"Marked: chew @ {current_ms/1000:.2f}s")

            elif key in UNDO_KEYS:
                if not history:
                    flash("No items can be deleted")
                else:
                    last = history.pop()
                    if last['type'] == 'chew':
                        r = last['row']
                        if len(rows[r]) > 1 and rows[r][-1] == last.get('ts'):
                            rows[r].pop()
                        elif len(rows[r]) > 1:
                            rows[r].pop()
                        flash("Deleted: last chew")
                    elif last['type'] == 'reach':
                        r = last['row']
                        if r == len(rows) - 1 and len(rows[r]) == 1:
                            rows.pop()
                            flash("Deleted: last new_chip")
                        else:
                            flash("Cannot delete (has data in row)")

            elif key == BACKWARD_KEY:
                cap.set(cv2.CAP_PROP_POS_MSEC, max(0, current_ms - 1000))
                ret, f2 = cap.read()
                if ret:
                    frame = f2
                flash("-1 Sec")

            elif key == FORWARD_KEY:
                cap.set(cv2.CAP_PROP_POS_MSEC, min(total_ms, current_ms + 1000))
                ret, f2 = cap.read()
                if ret:
                    frame = f2
                flash("+1 Sec")

            elif key == SMALL_BACK_KEY:
                cap.set(cv2.CAP_PROP_POS_MSEC, max(0, current_ms - 200))
                ret, f2 = cap.read()
                if ret:
                    frame = f2
                flash("-0.2 Sec")

            elif key == SMALL_FWD_KEY:
                cap.set(cv2.CAP_PROP_POS_MSEC, min(total_ms, current_ms + 200))
                ret, f2 = cap.read()
                if ret:
                    frame = f2
                flash("+0.2 Sec")

            elif key == FRAME_BACK_KEY:
                pos = cap.get(cv2.CAP_PROP_POS_FRAMES)
                cap.set(cv2.CAP_PROP_POS_FRAMES, max(0, pos - 2))
                ret, f2 = cap.read()
                if ret:
                    frame = f2

            elif key == FRAME_FWD_KEY:
                ret, f2 = cap.read()
                if ret:
                    frame = f2

            elif key in SPEED_KEYS:
                speed = SPEED_KEYS[key]
                flash(f"Speed x{speed:g}")

            elif key == PAUSE_KEY:
                paused = not paused

            elif key == NEXT_KEY:
                cap.release()
                return rows, False

            elif key == QUIT_KEY:
                cap.release()
                return rows, True

            # Autosave: write the partial result to CSV after every click (survives a forced quit)
            if marked and autosave_dir:
                try:
                    export_rows(video_path, rows, autosave_dir)
                except Exception:
                    pass

        if not paused and key == 255:
            ret, f2 = cap.read()
            if not ret:
                break
            frame = f2

    cap.release()
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
    ap = argparse.ArgumentParser(description='Video annotator (GUI buttons + slow playback version)')
    ap.add_argument('folder', nargs='?', default='.', help='video folder (default: current folder)')
    ap.add_argument('--outdir', default='csv_outputs', help='output folder (one CSV per video)')
    args = ap.parse_args()

    videos = find_videos(args.folder)
    if not videos:
        print(f"No videos found in: {args.folder}")
        sys.exit(1)

    for vp in videos:
        print(f"Playing: {vp}")
        rows, quit_all = [], False
        try:
            rows, quit_all = play_and_collect(vp, autosave_dir=args.outdir)
        except KeyboardInterrupt:
            print("(interrupted — the autosaved partial results remain in outdir)")
            quit_all = True
        out = export_rows(vp, rows, args.outdir)
        if out:
            print(f"Saved: {out}")
        else:
            print("(nothing new saved — if an autosaved file exists, it is the latest version)")
        if quit_all:
            break
    cv2.destroyAllWindows()


if __name__ == '__main__':
    main()
