#!/usr/bin/env bash
# Build the paper and report anything a referee would notice.
# Run from the research folder (MagicSquareOfSquaresResearch):
#   bash ../magic-square-tools/build.sh [name]   (default: magic_square_of_squares)
# Reads TEX/<name>.tex, writes PDF/<name>.pdf and build.log.
# PDF/MagicSquaresOfSquaresResearch.pdf is the original document that started
# the project and is never rebuilt or overwritten.
set -uo pipefail
NAME="${1:-magic_square_of_squares}"
PROTECTED="MagicSquaresOfSquaresResearch"
ROOT="$(pwd)"
LOG="$ROOT/build.log"

command -v latexmk >/dev/null || { echo "latexmk not found. Install TeX Live / MacTeX / MiKTeX."; exit 1; }
if [ "$NAME" = "$PROTECTED" ]; then
  echo "refusing to build $NAME: PDF/$PROTECTED.pdf is the original and stays unchanged"; exit 1
fi
[ -f "TEX/$NAME.tex" ] || { echo "TEX/$NAME.tex not found (run from the research folder)"; exit 1; }

cd TEX
rm -f "$NAME".{aux,log,out,fls,fdb_latexmk,toc} "$LOG"
latexmk -pdf -interaction=nonstopmode "$NAME.tex" > "$LOG" 2>&1
STATUS=$?

ERR=$(grep -c '^!' "$LOG")
OVER=$(grep -c 'Overfull' "$LOG")
UNDER=$(grep -c 'Underfull' "$LOG")
UNRES=$(awk '/Run number 2/,0' "$LOG" | grep -c 'undefined')

echo "build exit      : $STATUS"
echo "latex errors    : $ERR"
echo "overfull boxes  : $OVER"
echo "underfull boxes : $UNDER"
echo "unresolved refs : $UNRES   (after final pass)"

if [ "$ERR" -gt 0 ]; then
  echo; echo "--- errors ---"; grep -n -A3 '^!' "$LOG" | head -40
fi
if [ "$OVER" -gt 0 ]; then
  echo; echo "--- overfull locations ---"; grep 'Overfull' "$LOG" | sort -u | head -20
fi

if [ -f "$NAME.pdf" ]; then
  mkdir -p "$ROOT/PDF"
  mv -f "$NAME.pdf" "$ROOT/PDF/$NAME.pdf"
  echo "PDF written to PDF/$NAME.pdf ($(grep -o 'Output written on .*' "$LOG" | tail -1 | grep -o '([^)]*)'))"
fi
latexmk -c "$NAME.tex" >/dev/null 2>&1
[ "$ERR" -eq 0 ] && [ "$STATUS" -eq 0 ]
