# Reproducibility notes

## Scope

This repository reproduces the experimental results in Section 9 of the submitted manuscript. It is intentionally organized by paper result rather than by internal project phases.

## Determinism

The ordering-stress campaign uses fixed seeds `2026092100` through `2026092119`. Frozen random comparison orders are generated deterministically from seed offsets.

## Three independent numerical paths

- `certificate()` constructs the reduced A-only return operator.
- `sweep_operator()` independently assembles the full sweep propagator.
- `ibmi_sweep()` performs actual block inverse reconstruction.

This separation is a reproducibility guard against circular validation.

## Floating-point floor

For actual inverse reconstruction, the public-SPD runner uses

```text
max(1e-14, 100 * machine_epsilon * cond2(A))
```

and estimates local ratios only while the destination error is safely above the floor.

## Public data

SuiteSparse matrices are downloaded rather than redistributed. The frozen targets are `HB/bcsstk01` through `HB/bcsstk08`.

## Platform variation

Eigenvalues, condition estimates, timings, and last digits can vary with Python/NumPy/SciPy and BLAS/LAPACK builds. Use `scripts/capture_environment.py` to record the execution environment.

## Reviewer verification

`paper/EXPECTED_OUTPUTS.md` contains paper-level checkpoints. `scripts/verify_paper_outputs.py` checks the generated CSV files against robust tolerances rather than demanding bitwise equality.
