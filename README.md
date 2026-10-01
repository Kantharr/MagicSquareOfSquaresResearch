# Magic Square of Squares Tools

Companion computational toolchain for John Odom's research on the 3x3 magic
square of squares problem: nine distinct perfect squares in a 3x3 grid whose
rows, columns and both diagonals have the same sum.

This repo holds scripts only. The paper source (`TEX/`), the PDF (`PDF/`),
the research notes, the cited papers and all run output live in the research
folder `C:\Users\Owner\source\repos\MagicSquareOfSquaresResearch`.

## Layout

```
.
├── build.sh          -- build the paper's PDF (run from the research folder)
├── check.py          -- static consistency checks on the LaTeX source
├── verify_grid.py    -- check any 3x3 grid: line sums, squares, distinctness
└── center_search.py  -- search by centre E^2 via right triangles with hypotenuse E
```

## Requirements

- Python 3 (standard library only, so far).
- `latexmk` (MiKTeX or TeX Live) for `build.sh`.
- Later routes may need `gcc`, PARI/GP (`PARI_GP_PATH` if `gp` is not on
  `PATH`) and SageMath, as in the perfect-cuboid tools.

## Scripts

### `build.sh`

Run from the research folder. Reads `TEX/<name>.tex`, writes
`PDF/<name>.pdf` and `build.log`, and prints error, overfull/underfull-box
and unresolved-reference counts. Default name `MagicSquaresOfSquaresResearch`.

```sh
cd MagicSquareOfSquaresResearch
bash ../magic-square-tools/build.sh
```

### `check.py`

Same checker as the perfect-cuboid tools: dangling `\ref`s, duplicate labels,
uncited/undefined bibliography entries, unbalanced theorem environments,
statements without proofs, and an advisory on reused single-letter variables.

```sh
python ../magic-square-tools/check.py TEX/MagicSquaresOfSquaresResearch.tex
```

### `verify_grid.py`

With no arguments, checks two reference grids: Bremner's 1999 magic square
with 7 square entries (centre 425^2, magic sum 541875) and the Parker square
(9 squares, only 6 distinct, one diagonal fails). With nine arguments
(`373^2` style accepted), checks that grid. Exit 0 only for a genuine
solution.

### `center_search.py N`

Every 3x3 magic square has magic sum 3e (e the centre) and
a+i = b+h = c+g = d+f = 2e. With centre E^2, each of these pairs is
X^2, Y^2 with X^2 + Y^2 = 2E^2; setting X = p-q, Y = p+q turns this into
p^2 + q^2 = E^2 with step 2pq. So with

    T(E) = { p*q : p > q > 0, p^2 + q^2 = E^2 }

a solution with centre E^2 exists exactly when T(E) contains
s, r, r+s, r+2s for some r != s. The script builds T(E) for all E <= N from
Euclid's formula and reports primitive solutions and primitive 3-of-4 near
misses (7-square magic squares made of the centre and three full pairs).

```sh
python ../magic-square-tools/center_search.py 200000 --quiet   # about 5 s
```

Result so far: no solutions and no 3-of-4 near misses for E <= 200000.
Exit code is 0 only if a solution is found.
