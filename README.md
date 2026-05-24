# solar2d1d

MATLAB code for coupled 2D-optical / 1D-electrical simulation of thin-film
solar cells. The optical model is a rigorous coupled-wave analysis (RCWA)
solver for periodic structures; the electrical model is a drift-diffusion
Newton solver on a discontinuous-Galerkin / Lobatto polynomial basis. A
differential-evolution optimiser drives device-design studies.

Originally written as part of Tom Anderson's PhD work. Currently being
tidied, tested, and migrated to OOP.

## Layout

| Path | Contents |
|---|---|
| `Run*.m`, `OptoElec.m` | Top-level driver / optimisation entry points |
| `Design*.m` | Builds the `sim` configuration struct (junction, sim setup, DE) |
| `RCWA/` | Optical (2D, RCWA) solver |
| `DDNewton/` | Electrical (1D, drift-diffusion + Newton) solver |
| `DDNewton/PolyFuncs/` | Lobatto polynomial basis operations |
| `Materials/` | Material database + `Make*.m` builders for refractive-index data |
| `Materials/DataFiles/` | Source spreadsheets/CSVs for material properties |
| `Conversions/` | One-line unit conversions |
| `DEA/` | Differential-evolution optimisation (Buehren) |
| `Saved Inputs/` | Sample `Design*.m` configurations |
| `Slaves/`, `*/Slaves/` | Empty by default; file-semaphore mailbox for DE workers |
| `GUI.m`, `GUI.fig` | Legacy GUIDE GUI (GUIDE removed in MATLAB R2025b) |
| `legacy/scratch/` | Historical scratch / dev-note files preserved for reference |
| `legacy/notebooks/` | Mathematica derivation notebooks (reference only) |
| `ReadMe.tex` | Original LaTeX notes on the algorithm |

## Running

A single forward simulation:

```matlab
RunSim
```

Driver scripts live at the repo root. A parallel optimisation run can be
launched with `runview.sh`, which `nohup`s one optimisation master plus
three DE slaves.

## Tests

Unit tests live under `tests/` and use MATLAB's `matlab.unittest`
framework. From the repo root:

```matlab
runtests('tests')          % discover and run everything under tests/
runAllTests                % same, but adds the repo to the path first
```

CI runs the suite on every push and PR via
`.github/workflows/matlab-tests.yml` using
[`matlab-actions/setup-matlab`](https://github.com/matlab-actions/setup-matlab)
on a GitHub-hosted Linux runner. JUnit + Cobertura coverage reports are
uploaded as workflow artefacts. The workflow needs no licence token for
public repositories; for a private repo, set the `MATHWORKS_TOKEN`
secret on the repo.

Current coverage is intentionally narrow: the unit tests target leaf
functions (`Conversions/`, `DDNewton/PolyFuncs/`, RCWA `sortmat`) that
take primitive arguments. Coverage will broaden once the `sim` struct is
migrated to classes.

## Status

Active modernisation:

1. Hygiene pass: remove editor backups, transient slave mailbox files,
   per-host result/timing artefacts; add `.gitignore`; move scratch
   files to `legacy/`; lift the hard-coded `switch 3` bandgap-profile
   selector in `OptoElec.m` onto `params`.
2. Data purge: drop `Results/` historical outputs (84 files, never
   referenced from code); drop duplicate `SSAM15G.xlsx` at root; move
   Mathematica notebooks to `legacy/notebooks/`; move paper-plot
   script + `prometheus1.csv` to `legacy/scratch/`; relocate root
   `CZTS.xlsx` into `Materials/DataFiles/` and fix the read path.
3. Golden regression tests (planned).
4. Migrate the `sim` struct to MATLAB classes; replace GUIDE GUI with
   App Designer.

## Known issues

- `GUI.m`/`GUI.fig` rely on GUIDE, which is removed in MATLAB R2025b.
- `OptoElec.m` `catch` swallows all errors and writes NaNs to the
  figure-of-merit; real failures are silent.
- `RCWA.m` silently falls back to `pinv` when `rcond < 1e-10`; ill-
  conditioned slices are masked rather than reported.
- No automated tests yet.
