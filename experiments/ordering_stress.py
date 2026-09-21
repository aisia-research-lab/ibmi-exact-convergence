from pathlib import Path
import csv,itertools,sys
import numpy as np
import scipy.linalg as la
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from ibmi_cert import *
OUT=Path(__file__).resolve().parents[1]/"results"; OUT.mkdir(exist_ok=True)
SEEDS=list(range(2026092100,2026092120))

def precision(n,seed,regime):
    rng=np.random.default_rng(seed); W=np.zeros((n,n))
    if regime=="sparse":
        for i in range(n):
            for j in range(i+1,n):
                if rng.random()<.35: W[i,j]=W[j,i]=rng.uniform(.01,.16)
    elif regime=="dense":
        X=rng.uniform(.005,.08,(n,n)); W=np.triu(X,1); W=W+W.T
    elif regime=="heterogeneous":
        X=np.exp(rng.uniform(np.log(.002),np.log(.18),(n,n))); W=np.triu(X,1); W=W+W.T
    else: raise ValueError(regime)
    mx=np.max(np.sum(abs(W),axis=1))
    if mx>=.75: W*=.75/mx
    S=rng.choice([-1.,1.],(n,n)); S=np.triu(S,1); S=S+S.T
    H=np.eye(n)+W*S; la.cholesky(H,lower=True); return H

def interlaced_blocks(n,K=6):
    x=np.arange(n); return [np.setdiff1d(x,x[x%K==k]) for k in range(K)]

def main(quick=False):
    seeds=SEEDS[:2] if quick else SEEDS; rows=[]
    for regime in ["sparse","dense","heterogeneous"]:
        for seed in seeds:
            H=precision(8,seed,regime); O=list(itertools.permutations(range(8)))
            V=np.array([singleton_cycle_factor(H,o) for o in O]); native=singleton_cycle_factor(H,tuple(range(8)))
            rng=np.random.default_rng(seed+77); R=[tuple(rng.permutation(8)) for _ in range(20)]
            br=min(singleton_cycle_factor(H,o) for o in R); _,tw=two_opt_singleton(H,tuple(range(8))); best=float(V.min())
            rows.append(["singleton",regime,seed,8,native,best,float(np.median(V)),float(V.max()),br,tw,native/max(best,1e-300),br/max(best,1e-300),tw/max(best,1e-300)])
    for regime in ["sparse","dense","heterogeneous"]:
        for seed in seeds:
            H=precision(48,seed,regime); A=la.inv(H); B=interlaced_blocks(48,6); scores=exhaustive_fixed_block_scores(A,B)
            vals=np.array([v for _,v in scores]); native=dict(scores)[tuple(range(6))]; best=float(vals.min())
            rng=np.random.default_rng(seed+88); R=[tuple(rng.permutation(6)) for _ in range(20)]
            D=dict(scores); br=min(D[o] for o in R)
            rows.append(["general",regime,seed,6,native,best,float(np.median(vals)),float(vals.max()),br,np.nan,native/max(best,1e-300),br/max(best,1e-300),np.nan])
    with open(OUT/"ordering_stress.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["family","regime","seed","K","native","best","median","worst","best_random20","two_opt","native_best_ratio","random20_best_ratio","two_opt_best_ratio"]); w.writerows(rows)
    print("WROTE",OUT/"ordering_stress.csv",len(rows),"rows")
if __name__=="__main__": main("--quick" in sys.argv)
