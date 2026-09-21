from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; R=ROOT/"results"; F=ROOT/"figures"; F.mkdir(exist_ok=True)

p=R/"covariance_certificate.csv"
if p.exists():
    d=pd.read_csv(p).dropna(subset=["assembled_Q_factor"])
    fig,ax=plt.subplots(figsize=(4.4,4.0))
    x=np.log10(np.maximum(d.certificate,1e-300)); y=np.log10(np.maximum(d.assembled_Q_factor,1e-300))
    ax.scatter(x,y); lo=min(x.min(),y.min()); hi=max(x.max(),y.max()); ax.plot([lo,hi],[lo,hi])
    ax.set_xlabel(r"$\log_{10}\rho(M_A)^2$"); ax.set_ylabel(r"$\log_{10}\rho(Q)^2$"); ax.set_title("Exact certificate vs assembled sweep")
    fig.tight_layout(); fig.savefig(F/"certificate_vs_sweep.pdf"); plt.close(fig)

files=[R/f"trajectory_bcsstk03_overlap{x:02d}_native.csv" for x in [0,5,10,20]]
if all(p.exists() for p in files):
    fig,ax=plt.subplots(figsize=(5.0,3.8))
    for ov,p in zip([0,5,10,20],files):
        d=pd.read_csv(p); ax.semilogy(d.sweep,d.actual_rel_error2,label=f"{ov}% overlap")
    ax.set_xlabel("Sweep"); ax.set_ylabel("Relative inverse error (2-norm)"); ax.legend()
    fig.tight_layout(); fig.savefig(F/"bcsstk03_trajectories.pdf"); plt.close(fig)

p=R/"ordering_stress.csv"
if p.exists():
    d=pd.read_csv(p); d=d[d.family=="general"]
    vals=[np.log10(d[d.regime==r].native_best_ratio) for r in ["dense","heterogeneous","sparse"]]
    fig,ax=plt.subplots(figsize=(5.0,3.8)); ax.boxplot(vals,tick_labels=["dense","heterogeneous","sparse"])
    ax.set_ylabel(r"$\log_{10}(\gamma_{\rm native}/\gamma_{\rm best})$"); ax.set_xlabel("Frozen precision family")
    fig.tight_layout(); fig.savefig(F/"ordering_sensitivity.pdf"); plt.close(fig)
print("Figures written to",F)
