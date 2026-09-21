from pathlib import Path
import csv,itertools,sys
import numpy as np
import scipy.linalg as la
from scipy.io import mmread
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from ibmi_cert import *

ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/"data"; OUT=ROOT/"results"; OUT.mkdir(exist_ok=True)
EPS=np.finfo(float).eps; TOL=1e-8

def robust_ratio(e,floor):
    e=np.asarray(e); ratios=e[1:]/np.maximum(e[:-1],1e-300); good=np.where(e[1:]>10*floor)[0]
    if len(good)>=3:
        idx=good[-min(5,len(good)):]; return float(np.median(ratios[idx]))
    return np.nan

def trajectory(A,B,order,max_sweeps=300):
    H=la.inv(A); hnorm=la.norm(H,2); cond=float(np.linalg.cond(A))
    floor=max(1e-14,100*EPS*cond); X=np.eye(A.shape[0]); E=X-H; Q=sweep_operator(A,B,order)
    rows=[]; ae=[]; oe=[]
    for s in range(max_sweeps+1):
        a=la.norm(X-H,2)/hnorm; o=la.norm(E,2)/hnorm
        rows.append([s,a,o]); ae.append(a); oe.append(o)
        if a<TOL and o<TOL: break
        X=ibmi_sweep(A,X,B,order); E=Q@E@Q.T
    return rows,robust_ratio(oe,1e-14),robust_ratio(ae,floor),cond,floor

def main():
    summary=[]
    for i in range(1,9):
        name=f"bcsstk{i:02d}"; p=DATA/f"{name}.mtx"
        if not p.exists(): print("MISSING",p); continue
        A=mmread(p); A=A.toarray() if hasattr(A,"toarray") else np.asarray(A)
        A=scale_diagonal_congruence((A+A.T)/2)
        for ov in [0,.05,.10,.20]:
            B=contiguous_overlapping_blocks(A.shape[0],4,ov); C=cache_local_maps(A,B)
            orders=list(itertools.permutations(range(4)))
            best=min(orders,key=lambda o:certificate_from_cache(C,o))
            rng=np.random.default_rng(20260921+A.shape[0]+int(100*ov)); rnd=tuple(rng.permutation(4))
            for label,o in [("native",(0,1,2,3)),("random",rnd),("cgso",best)]:
                gamma=certificate_from_cache(C,o); tr,op,act,cond,floor=trajectory(A,B,o)
                fn=OUT/f"trajectory_{name}_overlap{int(100*ov):02d}_{label}.csv"
                with open(fn,"w",newline="") as f:
                    w=csv.writer(f); w.writerow(["sweep","actual_rel_error2","operator_rel_error2"]); w.writerows(tr)
                summary.append([name,ov,label,gamma,op,act,len(tr)-1,cond,floor])
    with open(OUT/"public_spd_summary.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["matrix","overlap","order","certificate","operator_pre_floor_ratio","actual_pre_floor_ratio","sweeps","cond2_scaled","floor_guard"]); w.writerows(summary)
    print("WROTE",OUT/"public_spd_summary.csv",len(summary),"rows")
if __name__=="__main__": main()
