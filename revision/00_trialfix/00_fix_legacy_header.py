#!/usr/bin/env python3
"""
00_fix_legacy_header.py  (Trial-column fix, step -1: header normalisation of the raw data; v1, 2026-09-15)

One annotation file, annotation_raw_data/solo_Human/20170706_07_03_B.csv, was exported by an older
version of the annotation tool and carries non-English column headers
    <U+56DE U+6578>, 1<U+565B U+307F U+76EE>, 2<...>, ..., 31<...>   (i.e. "trial number", "k-th bite")
whereas every other file uses
    Trial, Chew_1, Chew_2, ...
This script rewrites ONLY the header line of such files to the standard names
(trial-number column -> Trial, k-th-bite columns -> Chew_k). The data rows are copied byte for byte (line endings included), so the
event times are unchanged; 00_make_chewonly_data.py then treats the file like all the others.

  - The original file is copied to 00_trialfix/legacy_header_backup/<relative path> before it is rewritten
    (nothing is deleted).
  - A report (legacy_header_report.txt, next to this script) lists old/new header and the MD5 of the data
    rows before and after (must be identical).
  - Re-running is safe: files whose header is already standard are left untouched.

Usage:
  python 00_fix_legacy_header.py --src ../annotation_raw_data            # rewrite
  python 00_fix_legacy_header.py --src ../annotation_raw_data --dry_run  # only report
"""
import argparse
import hashlib
import os
import re
import shutil

LEGACY_TRIAL = "\u56de\u6578"            # legacy "trial number" header (2 CJK characters)
LEGACY_CHEW = re.compile(r"^(\d+)\u565b\u307f\u76ee$")   # legacy "k-th bite" header (digits + 3 CJK characters)


def map_header(name):
    n = name.strip().lstrip("\ufeff")
    if n == LEGACY_TRIAL:
        return "Trial"
    m = LEGACY_CHEW.match(n)
    if m:
        return f"Chew_{int(m.group(1))}"
    return None


def split_header(raw):
    for eol in (b"\r\n", b"\n"):
        i = raw.find(eol)
        if i >= 0:
            return raw[:i], eol, raw[i + len(eol):]
    return raw, b"", b""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="../annotation_raw_data")
    ap.add_argument("--dry_run", action="store_true")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    bk_root = os.path.join(here, "legacy_header_backup")
    lines = [f"legacy header normalisation  src={os.path.abspath(a.src)}  dry_run={a.dry_run}"]
    n_fixed = 0
    for r, ds, fs in os.walk(a.src):
        for f in sorted(fs):
            if not f.lower().endswith(".csv"):
                continue
            p = os.path.join(r, f)
            raw = open(p, "rb").read()
            head, eol, body = split_header(raw)
            try:
                cols = head.decode("utf-8").split(",")
            except UnicodeDecodeError:
                continue
            mapped = [map_header(c) for c in cols]
            if not any(mapped):
                continue                       # header already standard (or unrelated)
            new_cols = [m if m else c for c, m in zip(cols, mapped)]
            rel = os.path.relpath(p, a.src)
            md5_body = hashlib.md5(body).hexdigest()
            lines.append(f"[FIX] {rel}")
            lines.append(f"      old header: {','.join(cols)}")
            lines.append(f"      new header: {','.join(new_cols)}")
            lines.append(f"      data rows : {body.count(eol) if eol else 0} lines, md5 {md5_body} (unchanged)")
            if not a.dry_run:
                bk = os.path.join(bk_root, rel)
                os.makedirs(os.path.dirname(bk), exist_ok=True)
                if not os.path.exists(bk):
                    shutil.copy2(p, bk)
                new_raw = ",".join(new_cols).encode("utf-8") + eol + body
                open(p, "wb").write(new_raw)
                assert split_header(open(p, "rb").read())[2] == body
                lines.append(f"      original saved to 00_trialfix/legacy_header_backup/{rel}")
            n_fixed += 1
    lines.append(f"files with legacy headers: {n_fixed}" + ("" if n_fixed else " (nothing to do)"))
    rep = "\n".join(lines) + "\n"
    print(rep, end="")
    if not a.dry_run:
        open(os.path.join(here, "legacy_header_report.txt"), "w", encoding="utf-8").write(rep)


if __name__ == "__main__":
    main()
