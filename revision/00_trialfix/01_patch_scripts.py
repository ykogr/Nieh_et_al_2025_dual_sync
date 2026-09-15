#!/usr/bin/env python3
"""
01_patch_scripts.py  (Trial-column fix, step 1: code-level fix)

Purpose
  Scan the Python scripts in the given folders and fix, in place (after taking
  a backup), two defects inherited from the old 03_Micro_Analysis_NOLPM.py.

  [EVENT-COLUMN PATCH]  load_times_from_csv()
      Old: every CSV cell was numericised -> the Trial column (1, 2, ..., k) leaked
          into the event series as spurious events at "t = 1 s, 2 s, ...".
      New: use only the columns whose name starts with "chew" (long format with a "t" column is handled as before).

  [FLOAT-TOL PATCH]  NLOPM_match()
      Old: abs(a - b) <= delta was compared in raw floating point.
          For times on the 0.1 s grid, 0.3 - 0.2 = 0.09999999999999998 and
          0.7 - 0.6 = 0.09999999999999998 but 0.8 - 0.7 = 0.10000000000000009,
          so events exactly ±delta apart were matched or not depending on
          their values (with delta = 0.1 the ±0.1 s matches were unstable).
      New: compare with an added tolerance FLOAT_TOL = 1e-6 s (far below the 0.1 s grid).
          The same tolerance is applied to the look-ahead equal-distance test (<=),
          making the original rule "on a tie, prefer the later candidate" deterministic.

The target functions in 02_Micro_Analysis_NLOPM_periods*.py,
03_Micro_Analysis_NOLPM*.py, 01_null_symmetry_check*.py,
02_jitter_surrogates*.py, 03_reactive_shift_test*.py etc. were all copied from
the same prototype, so they can be replaced with one identical pattern.
A file that contains NLOPM_match / load_times_from_csv but where the pattern is
not found is reported as "manual check required" (nothing is changed).

Usage (run inside iscience_revision_analysis/00_trialfix/)
  python 01_patch_scripts.py ../gridnull_main ../R1_2 ../R2_2_2_5_2_10_3_2_3_3 ../R2_4_3_1
  python 01_patch_scripts.py --dry_run ../R1_2      # only show what would change

  The pre-change file is kept alongside as <original file name>.pre_trialfix.bak.
  Files that are already patched (marker present) are skipped (re-running is safe).

  [GRID-NULL PATCH]  R2-4_zS_reproduction.py only (v4.2)
      Old: circular-shift offsets drawn from a continuous uniform U(0,T) (the old
          convention acknowledged as an error in the letter's General note).
      New: adds --null {continuous,grid}. With grid, shifts are integer multiples of
          the 0.1 s grid only (same construction as R2-4_zS_nullcheck.py). The default
          stays continuous; run_R2_4_trialfix.sh passes --null grid.

  [OBSERVED-ARG PATCH]  R3-1_power_sensitivity.R only (v4.2)
      Old: the closing "observed P3 contrast = ..." line of the output was a hard-coded string.
      New: prints the values of the --observed / --observed_se arguments (passed by kick from lmm_results.txt).
"""
import argparse
import os
import re
import shutil
import sys

MARK_EVENT = "[EVENT-COLUMN PATCH"
MARK_FLOAT = "[FLOAT-TOL PATCH"

NEW_LOADER = '''def load_times_from_csv(path: str) -> np.ndarray:
    """
    [EVENT-COLUMN PATCH 2026-09-07]
    Read chewing-event times from one annotation CSV (one participant).
      - Long format with a "t" column: use that column (unchanged behaviour).
      - Wide format (Trial, Chew_1, Chew_2, ...): use ONLY the columns whose
        name starts with "chew". The previous implementation numericised
        EVERY cell (df.stack / to_numpy().ravel()), so the Trial index
        (1, 2, ..., k) leaked into the event series as spurious events at
        t = 1 s, 2 s, ..., k s.
    Keep non-negative values only; sort and drop duplicates (as before).
    """
    df = _read_csv_smart(path)
    if "t" in df.columns:
        s = pd.to_numeric(df["t"], errors="coerce")
    else:
        cols = [c for c in df.columns if str(c).strip().lower().startswith("chew")]
        if not cols:
            # fallback: drop obvious index-like columns, keep the rest
            drop = {"trial", "bite", "index", "unnamed: 0", ""}
            cols = [c for c in df.columns if str(c).strip().lower() not in drop]
        s = pd.to_numeric(pd.Series(df[cols].to_numpy().ravel()), errors="coerce")
    arr = s.dropna().astype(float).values
    arr = arr[arr >= 0]
    arr = np.unique(np.sort(arr))
    return arr
'''

FLOAT_TOL_BLOCK = '''# [FLOAT-TOL PATCH 2026-09-07] tolerance for comparisons against delta.
# Event times are on a 0.1-s grid; in binary floating point the difference of
# two grid values is not exactly a multiple of 0.1 (e.g. 0.8 - 0.7 =
# 0.10000000000000009 > 0.1), so events exactly delta apart were matched or
# not depending on their float representation. 1e-6 s is far below the grid.
FLOAT_TOL = float(os.environ.get("NLOPM_FLOAT_TOL", "1e-6"))


'''

# (regex, replacement)  -- applied inside the NLOPM_match function body only
NLOPM_RULES = [
    (re.compile(r"(while j < n and b\[j\] < a\[i\] - delta)(:)"),
     r"\1 - FLOAT_TOL\2  # [FLOAT-TOL PATCH]"),
    (re.compile(r"(while i < m and a\[i\] < b\[j\] - delta)(:)"),
     r"\1 - FLOAT_TOL\2  # [FLOAT-TOL PATCH]"),
    (re.compile(r"\(b\[j_star \+ 1\] <= a\[i\] \+ delta\)"),
     r"(b[j_star + 1] <= a[i] + delta + FLOAT_TOL)"),
    (re.compile(r"\(abs\(a\[i\] - b\[j_star \+ 1\]\) <= abs\(a\[i\] - b\[j_star\]\)\)"),
     r"(abs(a[i] - b[j_star + 1]) <= abs(a[i] - b[j_star]) + FLOAT_TOL)"),
    (re.compile(r"(if abs\(a\[i\] - b\[j_star\]\) <= delta)(:)"),
     r"\1 + FLOAT_TOL\2  # [FLOAT-TOL PATCH]"),
]


# R2-4 reimplementation family (R2-4_zS_reproduction.py, R2-4_zS_nullcheck.py,
# R2-4_jitter_simulation.py): function `nlopm_matches`. Reads Chew columns only
# (no Trial leak), but the lower bound `B[j] < a - delta` and the tie tests have
# no tolerance, so matches exactly delta apart depend on float representation.
NLOPM2_RULES = [
    (re.compile(r"(while j < m and B\[j\] < a - delta)(:)"),
     r"\1 - FLOAT_TOL\2  # [FLOAT-TOL PATCH]"),
    (re.compile(r"while k < m and B\[k\] <= a \+ delta \+ 1e-12:"),
     r"while k < m and B[k] <= a + delta + FLOAT_TOL:  # [FLOAT-TOL PATCH]"),
    (re.compile(r"best = min\(cands, key=lambda kk: \(abs\(B\[kk\] - a\), B\[kk\]\)\)"),
     r"best = min(cands, key=lambda kk: (round(abs(B[kk] - a), 6), B[kk]))  # [FLOAT-TOL PATCH]"),
    (re.compile(r"if i \+ 1 < n and abs\(A\[i \+ 1\] - B\[best\]\) < abs\(a - B\[best\]\):"),
     r"if i + 1 < n and abs(A[i + 1] - B[best]) < abs(a - B[best]) - FLOAT_TOL:  # [FLOAT-TOL PATCH]"),
    (re.compile(r"best2 = min\(earlier, key=lambda kk: \(abs\(B\[kk\] - a\), B\[kk\]\)\)"),
     r"best2 = min(earlier, key=lambda kk: (round(abs(B[kk] - a), 6), B[kk]))  # [FLOAT-TOL PATCH]"),
]


# Sliding-window time axis (03_Micro_Analysis_NOLPM_gridnull*.py).
# Window centres were start + 2.5 + k. With the Trial artefact, `start` was
# 1.0 s for EVERY pair, so all pairs happened to share one grid (3.5, 4.5, ...).
# With real chew onsets (3.2-19.5 s, 0.1-s grid) the per-pair grids no longer
# coincide and the wide table / Bayesian bins / change-point group curve break.
# Default here: absolute video time, with the first window edge snapped up to
# the next whole second so every pair sits on the legacy x.5-s grid
# (3.5, 4.5, ... in video time; a pair starting at 4.3 s gets 7.5, 8.5, ...).
# Time since the first common chew is available via NLOPM_TIME_AXIS=relative.
MARK_TIME = "[TIME-AXIS PATCH"
TIME_AXIS_BLOCK = '''# [TIME-AXIS PATCH 2026-09-07] "absolute" (default): t_center is video time,
# window edges snapped to whole seconds so that all pairs share the legacy
# x.5-s grid; "relative": t_center = seconds since the first common chew
# (overlap start). Either way all pairs share one 1-s grid (needed downstream).
TIME_AXIS = os.environ.get("NLOPM_TIME_AXIS", "absolute")


'''
TIME_RULES = [
    (re.compile(r'^([ \t]*)centers = np\.arange\(start \+ win/2\.0, end - win/2\.0 \+ 1e-9, step\)\s*$', re.M),
     r'\1if TIME_AXIS == "absolute":  # [TIME-AXIS PATCH]' '\n'
     r'\1    _first = float(np.ceil(start - 1e-9))' '\n'
     r'\1    centers = np.arange(_first + win/2.0, end - win/2.0 + 1e-9, step)' '\n'
     r'\1else:' '\n'
     r'\1    centers = np.arange(start + win/2.0, end - win/2.0 + 1e-9, step)'),
    (re.compile(r'^([ \t]*)obs_df = pd\.DataFrame\(obs_rows, columns=\["t_center", "S", "nA", "nB", "nM"\]\)\s*$', re.M),
     r'\1obs_df = pd.DataFrame(obs_rows, columns=["t_center", "S", "nA", "nB", "nM"])' '\n'
     r'\1# [TIME-AXIS PATCH] keep video time in t_center_abs; t_center is on a shared 1-s grid' '\n'
     r'\1obs_df["t_center_abs"] = np.round(obs_df["t_center"], 1)' '\n'
     r'\1if TIME_AXIS == "absolute":' '\n'
     r'\1    obs_df["t_center"] = obs_df["t_center_abs"]' '\n'
     r'\1else:' '\n'
     r'\1    obs_df["t_center"] = np.round(obs_df["t_center"] - start, 1)'),
    (re.compile(r'"t_center": row\["t_center"\], "S": row\["S"\],'),
     r'"t_center": row["t_center"], "t_center_abs": row["t_center_abs"], "S": row["S"],'),
]
TIME_RULE_EXPECT = [1, 1, 2]   # expected match counts per rule inside sliding_window_with_dt


MARK_PERIOD = "[PERIOD-ORIGIN PATCH"
PERIOD_BLOCK = '''# [PERIOD-ORIGIN PATCH 2026-09-08] "absolute" (default): periods t0-t1 are video
# time, intersected with the pair's overlap interval; "relative": seconds since
# the first shared chew (overlap start), which is what the code did before.
PERIOD_ORIGIN = os.environ.get("NLOPM_PERIOD_ORIGIN", "absolute")
'''
PERIOD_RX = re.compile(
    r"^([ \t]*)seg_start = start \+ t0\n[ \t]*seg_end = min\(start \+ t1, end\)\n", re.M)
PERIOD_REPL = (
    r"\1if PERIOD_ORIGIN == 'absolute':  # [PERIOD-ORIGIN PATCH]\n"
    r"\1    seg_start = max(float(t0), start)\n"
    r"\1    seg_end = min(float(t1), end)\n"
    r"\1else:\n"
    r"\1    seg_start = start + t0\n"
    r"\1    seg_end = min(start + t1, end)\n")



# ---------------------------------------------------------------------------
# [GRID-NULL PATCH]  R2-4_zS_reproduction.py only.
# The reproduction script drew circular-shift offsets from a continuous
# uniform distribution (u ~ U(0, T)); the letter's General note replaces that
# convention by grid-aligned shifts (integer multiples of the 0.1-s grid).
# Adds  --null {continuous,grid}  (default: continuous = original behaviour;
# run_R2_4_trialfix.sh passes --null grid) using the same shift construction
# as R2-4_zS_nullcheck.py, and records the setting in the summary header.
MARK_GRIDNULL = "[GRID-NULL PATCH"
GRIDNULL_RULES = [
    ("def zs_period(t_left, t_right, p0, p1, delta, n_perm, rng):",
     "def zs_period(t_left, t_right, p0, p1, delta, n_perm, rng,\n"
     "              null_shift=\"continuous\", grid=0.1):  # [GRID-NULL PATCH]"),
    ("    s_null = np.empty(n_perm)\n"
     "    for r in range(n_perm):\n"
     "        u = rng.uniform(0.0, T)\n"
     "        Bs = np.sort(tmin + (B - tmin + u) % T)\n",
     "    s_null = np.empty(n_perm)\n"
     "    n_steps = max(int(round(T / grid)), 2)  # [GRID-NULL PATCH]\n"
     "    for r in range(n_perm):\n"
     "        if null_shift == \"grid\":  # [GRID-NULL PATCH] same construction as R2-4_zS_nullcheck.py\n"
     "            u = rng.integers(1, n_steps) * grid\n"
     "            Bs = tmin + (B - tmin + u) % T\n"
     "            Bs = np.round(Bs / grid) * grid   # clean float error only\n"
     "            Bs = np.unique(Bs)\n"
     "        else:\n"
     "            u = rng.uniform(0.0, T)\n"
     "            Bs = np.sort(tmin + (B - tmin + u) % T)\n"),
    ("    ap.add_argument(\"--seed\", type=int, default=20260827)\n",
     "    ap.add_argument(\"--seed\", type=int, default=20260827)\n"
     "    ap.add_argument(\"--null\", choices=[\"continuous\", \"grid\"], default=\"continuous\",\n"
     "                    help=\"circular-shift offsets: continuous U(0,T) [original] or \"\n"
     "                         \"integer multiples of the grid [letter, General note]\")  # [GRID-NULL PATCH]\n"),
    ("                res = zs_period(tl, tr, p0, p1, args.delta,\n"
     "                                args.n_perm, rng)\n",
     "                res = zs_period(tl, tr, p0, p1, args.delta,\n"
     "                                args.n_perm, rng,\n"
     "                                null_shift=args.null, grid=args.grid)  # [GRID-NULL PATCH]\n"),
    ("        fh.write(f\"Permutations per cell : {args.n_perm}\\n\")\n",
     "        fh.write(f\"Permutations per cell : {args.n_perm}\\n\")\n"
     "        fh.write(f\"Null shift            : {args.null}\"\n"
     "                 + (f\" (offsets = k * {args.grid} s, grid-aligned)\" if args.null == \"grid\"\n"
     "                    else \" (offsets ~ U(0,T), original convention)\") + \"\\n\")  # [GRID-NULL PATCH]\n"),
]

# ---------------------------------------------------------------------------
# [OBSERVED-ARG PATCH]  R3-1_power_sensitivity.R only.
# The closing note of the power summary ("observed P3 contrast in the
# manuscript = ...") was a hard-coded string, not a computed value. Replace it
# by  --observed <est> --observed_se <se>  arguments (written by
# run_R2_4_trialfix.sh from ../R1_2/lmm_grid_out/lmm_results.txt).
MARK_OBSERVED = "[OBSERVED-ARG PATCH"
OBSERVED_NOTE_RX = re.compile(
    r'^cat\("\\nNote: observed P3 contrast in the (?:revised )?manuscript = [^"]*"\)\s*$', re.M)
OBSERVED_NOTE_NEW = (
    '# [OBSERVED-ARG PATCH] value comes from the command line, not from a hard-coded string\n'
    'observed    <- suppressWarnings(as.numeric(get_arg("--observed", NA)))\n'
    'observed_se <- suppressWarnings(as.numeric(get_arg("--observed_se", NA)))\n'
    'if (is.finite(observed)) {\n'
    '  cat(sprintf("\\nNote: observed P3 contrast (visible - invisible) in the re-run LMM = %+.3f%s.\\n",\n'
    '              observed, if (is.finite(observed_se)) sprintf(" (SE %.3f)", observed_se) else ""))\n'
    '  cat("      (from ../R1_2/lmm_grid_out/lmm_results.txt via R3-1_power_params.txt)\\n")\n'
    '} else {\n'
    '  cat("\\nNote: observed P3 contrast not supplied (--observed); see R3-1_power_params.txt.\\n")\n'
    '}')


def patch_gridnull(src):
    """R2-4_zS_reproduction.py: add --null grid. Returns (text, notes, manual)."""
    notes = []
    if MARK_GRIDNULL in src:
        return src, ["grid-null: already patched (skip)"], False
    out = src
    for old, new in GRIDNULL_RULES:
        if out.count(old) != 1:
            return src, [f"grid-null: pattern not found exactly once -> manual check: {old.strip().splitlines()[0][:60]}"], True
        out = out.replace(old, new)
    notes.append(f"grid-null: {len(GRIDNULL_RULES)} edits applied (--null {{continuous,grid}} added)")
    return out, notes, False


def patch_observed_arg(src):
    """R3-1_power_sensitivity.R: replace hard-coded note by --observed args."""
    if MARK_OBSERVED in src:
        return src, ["observed-arg: already patched (skip)"], False
    if len(OBSERVED_NOTE_RX.findall(src)) != 1 or "get_arg <- function" not in src:
        return src, ["observed-arg: hard-coded note not found exactly once -> manual check"], True
    out = OBSERVED_NOTE_RX.sub(lambda m: OBSERVED_NOTE_NEW, src)
    return out, ["observed-arg: hard-coded note replaced by --observed/--observed_se"], False


def find_function_span(lines, name):
    """Return (start, end) line indices [start, end) of top-level `def name(`."""
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith(f"def {name}("):
            start = i
            break
    if start is None:
        return None
    end = len(lines)
    for j in range(start + 1, len(lines)):
        ln = lines[j]
        if ln.strip() == "" or ln.startswith((" ", "\t")):
            continue
        end = j
        break
    # trim trailing blank lines back into the gap
    while end > start + 1 and lines[end - 1].strip() == "":
        end -= 1
    return start, end


def patch_text(src, fname):
    """Return (new_text, notes, needs_manual)."""
    notes = []
    manual = False
    lines = src.split("\n")

    # ---- loader ---------------------------------------------------------
    if "def load_times_from_csv(" in src:
        if MARK_EVENT in src:
            notes.append("loader: already patched (skip)")
        else:
            span = find_function_span(lines, "load_times_from_csv")
            if span is None or "_read_csv_smart" not in src:
                notes.append("loader: FOUND but unexpected form -> manual check")
                manual = True
            else:
                s, e = span
                lines = lines[:s] + NEW_LOADER.rstrip("\n").split("\n") + lines[e:]
                notes.append(f"loader: replaced lines {s + 1}-{e}")
    src2 = "\n".join(lines)

    # ---- NLOPM_match ----------------------------------------------------
    if "def NLOPM_match(" in src2:
        if MARK_FLOAT in src2:
            notes.append("NLOPM_match: already patched (skip)")
        else:
            lines = src2.split("\n")
            span = find_function_span(lines, "NLOPM_match")
            s, e = span
            body = "\n".join(lines[s:e])
            hits = 0
            for rx, rep in NLOPM_RULES:
                body, n = rx.subn(rep, body)
                hits += n
            if hits != len(NLOPM_RULES):
                notes.append(f"NLOPM_match: only {hits}/{len(NLOPM_RULES)} patterns matched -> manual check")
                manual = True
            else:
                new_body = body.split("\n")
                lines = lines[:s] + FLOAT_TOL_BLOCK.split("\n") + new_body + lines[e:]
                notes.append(f"NLOPM_match: {hits} comparisons patched (function at line {s + 1})")
                src2 = "\n".join(lines)
                if not re.search(r"^import os\b", src2, flags=re.M):
                    src2 = re.sub(r"^(import argparse\n)", r"\1import os\n", src2, count=1, flags=re.M)
                    if not re.search(r"^import os\b", src2, flags=re.M):
                        src2 = "import os\n" + src2
                    notes.append("added 'import os' (for NLOPM_FLOAT_TOL env override)")
    # ---- nlopm_matches (R2-4 reimplementation) ---------------------------
    if "def nlopm_matches(" in src2:
        if MARK_FLOAT in src2:
            notes.append("nlopm_matches: already patched (skip)")
        else:
            lines = src2.split("\n")
            s, e = find_function_span(lines, "nlopm_matches")
            body = "\n".join(lines[s:e])
            hits = 0
            for rx, rep in NLOPM2_RULES:
                body, n = rx.subn(rep, body)
                hits += n
            if hits != len(NLOPM2_RULES):
                notes.append(f"nlopm_matches: only {hits}/{len(NLOPM2_RULES)} patterns matched -> manual check")
                manual = True
            else:
                lines = lines[:s] + FLOAT_TOL_BLOCK.split("\n") + body.split("\n") + lines[e:]
                notes.append(f"nlopm_matches: {hits} comparisons patched (function at line {s + 1})")
                src2 = "\n".join(lines)
                if not re.search(r"^import os\b", src2, flags=re.M):
                    src2 = "import os\n" + src2
                    notes.append("added 'import os'")
    # ---- sliding-window time axis ------------------------------------------
    if "def sliding_window_with_dt(" in src2:
        if MARK_TIME in src2:
            notes.append("sliding_window_with_dt: already patched (skip)")
        else:
            lines = src2.split("\n")
            s, e = find_function_span(lines, "sliding_window_with_dt")
            body = "\n".join(lines[s:e])
            counts = []
            for rx, rep in TIME_RULES:
                body, n = rx.subn(rep, body)
                counts.append(n)
            if counts != TIME_RULE_EXPECT:
                notes.append(f"sliding_window_with_dt: pattern counts {counts} != {TIME_RULE_EXPECT} -> manual check")
                manual = True
            else:
                lines = lines[:s] + TIME_AXIS_BLOCK.split("\n") + body.split("\n") + lines[e:]
                notes.append(f"sliding_window_with_dt: time axis patched (function at line {s + 1})")
                src2 = "\n".join(lines)
                if not re.search(r"^import os\b", src2, flags=re.M):
                    src2 = "import os\n" + src2
                    notes.append("added 'import os'")
    # ---- period origin (R1_2 02_..._periods*.py / 03_pseudo_pairs*.py) -----
    # Periods were placed relative to the overlap start (seg_start = start + t0).
    # With the Trial artefact the overlap start was 1.0 s for every pair, so this
    # was effectively video time. Default now: video time (absolute), i.e. the
    # period [t0, t1) is intersected with the pair's overlap; NLOPM_PERIOD_ORIGIN=
    # relative restores "seconds since the first shared chew".
    if PERIOD_RX.search(src2):
        if MARK_PERIOD in src2:
            notes.append("period origin: already patched (skip)")
        else:
            src2, n = PERIOD_RX.subn(PERIOD_REPL, src2)
            head = src2.find("\ndef ")
            src2 = src2[:head] + "\n" + PERIOD_BLOCK + src2[head:]
            notes.append(f"period origin: {n} location(s) patched (video time by default)")
            if not re.search(r"^import os\b", src2, flags=re.M):
                src2 = "import os\n" + src2
                notes.append("added 'import os'")
    return src2, notes, manual


KIT_VERSION = "4.3.2"   # printed at start so that logs show which kit was used


def main():
    print(f"01_patch_scripts.py  (trialfix_kit v{KIT_VERSION}; patches: EVENT-COLUMN, FLOAT-TOL, TIME-AXIS, PERIOD-ORIGIN, GRID-NULL, OBSERVED-ARG)")
    ap = argparse.ArgumentParser()
    ap.add_argument("folders", nargs="+")
    ap.add_argument("--dry_run", action="store_true")
    ap.add_argument("--no_backup", action="store_true")
    args = ap.parse_args()

    total_changed = 0
    manual_files = []
    for folder in args.folders:
        folder = os.path.abspath(folder)
        if not os.path.isdir(folder):
            print(f"[WARN] not a folder: {folder}")
            continue
        print(f"\n=== {folder}")
        for fn in sorted(os.listdir(folder)):
            if not (fn.endswith(".py") or fn == "R3-1_power_sensitivity.R"):
                continue
            p = os.path.join(folder, fn)
            with open(p, encoding="utf-8", errors="replace") as f:
                src = f.read()
            if fn == "R3-1_power_sensitivity.R":
                new, notes, manual = patch_observed_arg(src)
            else:
                if ("def load_times_from_csv(" not in src and "def NLOPM_match(" not in src
                        and "def nlopm_matches(" not in src and "seg_start = start + t0" not in src):
                    continue
                new, notes, manual = patch_text(src, fn)
                if fn == "R2-4_zS_reproduction.py":
                    new, notes2, manual2 = patch_gridnull(new)
                    notes += notes2
                    manual = manual or manual2
            flag = "MANUAL" if manual else ("CHANGED" if new != src else "ok")
            print(f"  [{flag:7s}] {fn}")
            for n in notes:
                print(f"            - {n}")
            if manual:
                manual_files.append(p)
            if new != src and not args.dry_run:
                if not args.no_backup:
                    bak = p + ".pre_trialfix.bak"
                    if not os.path.exists(bak):
                        shutil.copy2(p, bak)
                with open(p, "w", encoding="utf-8") as f:
                    f.write(new)
                total_changed += 1
    print(f"\n[DONE] files changed: {total_changed}{' (dry run: nothing written)' if args.dry_run else ''}")
    if manual_files:
        print("[ATTENTION] the following files contain the target functions in an "
              "unexpected form and were NOT changed; please send them for manual patching:")
        for p in manual_files:
            print("   ", p)
        sys.exit(2)


if __name__ == "__main__":
    main()
