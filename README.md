# ibmi-exact-convergence

Reproducibility repository for:

**Binh T. Nguyen, _Exact Convergence Certificates and Ordering Theory for Overlapping Iterative Block Matrix Inversion_ (2026).**

**Manuscript status:** submitted to the **SIAM Journal on Matrix Analysis and Applications (SIMAX)**.  
**License:** MIT.

This repository is organized around the paper—not around the internal development history of the project. A reviewer should be able to clone the repository, run a smoke test, reproduce the synthetic experiments, download the eight public SuiteSparse SPD matrices, reproduce the public-matrix trajectories, and regenerate the three figures and Table 9.1 checkpoints.

## What this repository reproduces

| Paper result | Code | Main output |
|---|---|---|
| Theorems 5.1–5.2: exact A-only local transfers and return certificate | `src/ibmi_cert/core.py` | `certificate()` |
| Theorem 7.1: singleton precision-cycle factor | `src/ibmi_cert/ordering.py` | `singleton_cycle_factor()` |
| Proposition 7.4: factor-once exact permutation scoring | `src/ibmi_cert/ordering.py` | cached local maps |
| Section 9.1 / Figure 9.1 | `experiments/covariance_validation.py` | `results/covariance_certificate.csv` |
| Section 9.2 / Table 9.1 / Figure 9.2 | `experiments/public_spd_validation.py` | public-SPD summary + trajectories |
| Section 9.3 / Figure 9.3 | `experiments/ordering_stress.py` | `results/ordering_stress.csv` |
| Paper figures | `scripts/make_figures.py` | `figures/*.pdf` |

The implementation deliberately keeps three validation paths separate: the reduced exact certificate, an independently assembled full sweep operator, and actual IBMI block reconstruction. This makes the numerical validation non-circular.

## Quick start

```bash
git clone https://github.com/REPLACE-WITH-GITHUB-USERNAME/ibmi-exact-convergence.git
cd ibmi-exact-convergence

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt

python run_reproduction.py --mode smoke
```

Expected smoke result: all unit tests pass, a small covariance check runs, a small ordering-stress check runs, and figures that have sufficient inputs are regenerated.

## Full reviewer workflow

### A. Synthetic experiments: Figure 9.1 and Figure 9.3

```bash
python experiments/covariance_validation.py
python experiments/ordering_stress.py
python scripts/make_figures.py
python scripts/verify_paper_outputs.py
```

The full covariance campaign contains 30 cases: five covariance families (`EXP`, `RBF`, `IQUAD`, `M32`, `M52`), 1D and 2D geometries, and three sizes per geometry.

The frozen ordering-stress campaign uses seeds `2026092100`–`2026092119`. The general fixed-block experiment has `n=48`, `K=6`, and exhaustively scores all `6! = 720` permutations. The singleton-complement experiment has `n=8` and exhaustively scores all `8! = 40320` orders.

### B. Public SPD matrices: Table 9.1 and Figure 9.2

```bash
python scripts/download_suitesparse.py
python experiments/public_spd_validation.py
python scripts/make_figures.py
python scripts/verify_paper_outputs.py
```

This downloads `HB/bcsstk01` through `HB/bcsstk08` from the SuiteSparse Matrix Collection. The experiment applies diagonal congruence scaling and uses four contiguous blocks with nominal overlap `0%`, `5%`, `10%`, and `20%`. For every matrix/overlap pair it evaluates native, one frozen random, and certificate-minimizing orders.

The public matrices are downloaded rather than redistributed.

## Expected paper checkpoints

See [`paper/EXPECTED_OUTPUTS.md`](paper/EXPECTED_OUTPUTS.md) for reviewer-readable numerical checkpoints and [`paper/expected_table_9_1.csv`](paper/expected_table_9_1.csv) for the submitted Table 9.1 values.

In particular:

- **Figure 9.1:** certificate points should lie on the equality diagonal against the independently assembled sweep factor.
- **Table 9.1:** for `bcsstk03`, the native-order certificates at 0/5/10/20% overlap are approximately `0.983275 / 0.891011 / 0.678239 / 0.048749`.
- **Figure 9.2:** `bcsstk03` at 0% overlap does not reach tolerance within the 300-sweep cap; the 5/10/20% cases require 150/45/7 sweeps in the submitted run.
- **Figure 9.3:** the submitted median native-to-best factors for the general fixed-block stress campaign are approximately `6.76 / 19.66 / 31.09` for dense / heterogeneous / sparse precision families.

## Capture the execution environment

```bash
python scripts/capture_environment.py
```

This writes `results/environment.json` containing Python, OS/architecture, NumPy, SciPy, pandas, Matplotlib, CPU-count, and memory information.

The submitted experiments were CPU-only. Exact timings and final floating-point digits can vary with hardware and BLAS/LAPACK implementation.

## Numerical-floor guardrail

Actual inverse reconstruction can reach a conditioning-dependent floating-point floor before the exact asymptotic factor is visible. The public-SPD experiment therefore separates:

1. `gamma_A = rho(M_A)^2`, the exact spectral certificate;
2. independently propagated `E_{s+1} = Q E_s Q^T`;
3. actual IBMI reconstruction.

Rate estimation for actual reconstruction is restricted to errors safely above
`max(1e-14, 100 * eps * cond2(A))`.

## Repository layout

```text
.
├── CITATION.cff
├── LICENSE
├── Makefile
├── README.md
├── REPRODUCIBILITY.md
├── paper/
│   ├── manuscript.pdf
│   ├── EXPECTED_OUTPUTS.md
│   └── expected_table_9_1.csv
├── src/ibmi_cert/
│   ├── core.py
│   ├── ordering.py
│   └── problems.py
├── experiments/
│   ├── covariance_validation.py
│   ├── public_spd_validation.py
│   └── ordering_stress.py
├── scripts/
│   ├── capture_environment.py
│   ├── download_suitesparse.py
│   ├── make_figures.py
│   └── verify_paper_outputs.py
├── tests/
│   └── test_core.py
├── data/
├── results/
└── figures/
```

## License

Released under the [MIT License](LICENSE).

