from pathlib import Path
import pandas as pd, numpy as np, sys
ROOT=Path(__file__).resolve().parents[1]; R=ROOT/"results"
ok=True
def check(cond,msg):
    global ok
    print(("PASS " if cond else "FAIL ")+msg); ok &= bool(cond)

p=R/"covariance_certificate.csv"
if p.exists():
    d=pd.read_csv(p).dropna(subset=["assembled_Q_factor"])
    err=np.max(np.abs(d["log10_difference"]))
    check(err < 1e-6, f"Fig. 9.1 certificate/sweep agreement; max |log10 diff|={err:.3e}")
else: check(False,"missing covariance_certificate.csv")

p=R/"ordering_stress.csv"
if p.exists():
    d=pd.read_csv(p)
    g=d[d.family=="general"]
    check(len(g) in (6,60), f"ordering-stress general rows={len(g)}")
    if len(g)==60:
        for r,target in [("dense",6.76),("heterogeneous",19.66),("sparse",31.09)]:
            med=float(g[g.regime==r].native_best_ratio.median())
            check(abs(med-target)/target < .03, f"Fig. 9.3 {r} median={med:.4g}, paper≈{target}")
else: check(False,"missing ordering_stress.csv")

p=R/"public_spd_summary.csv"
if p.exists():
    d=pd.read_csv(p)
    n=d[d["order"]=="native"]
    check(len(n)==32, f"Table 9.1 native public-SPD rows={len(n)}")
    if len(n)==32:
        b=n[n.matrix=="bcsstk03"].sort_values("overlap")
        vals=b.certificate.to_numpy()
        exp=np.array([.983275,.891011,.678239,.048749])
        check(np.max(np.abs(vals-exp))<5e-6, "Table 9.1 bcsstk03 certificate checkpoints")
else:
    print("SKIP public SPD verification (SuiteSparse run not present)")
if not ok: sys.exit(1)
