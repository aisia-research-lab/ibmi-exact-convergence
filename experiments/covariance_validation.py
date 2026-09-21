from pathlib import Path
import csv, sys
import numpy as np
import scipy.linalg as la
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from ibmi_cert import *

OUT=Path(__file__).resolve().parents[1]/"results"; OUT.mkdir(exist_ok=True)

def main(quick=False):
    rows=[]
    kinds=["EXP","RBF","IQUAD","M32","M52"]
    sizes1=[64] if quick else [64,128,256]
    sizes2=[64] if quick else [64,144,256]
    for dim,sizes in [(1,sizes1),(2,sizes2)]:
        for n in sizes:
            P=points_1d(n) if dim==1 else points_2d(n)
            for kind in kinds:
                A=covariance_matrix(P,kind)
                B=contiguous_overlapping_blocks(n,4,.05)
                gamma=certificate(A,B)
                q2=np.nan
                if n<=144:
                    Q=sweep_operator(A,B); q2=float(max(abs(la.eigvals(Q)))**2)
                rows.append([kind,dim,n,gamma,q2,np.log10(max(gamma,1e-300))-np.log10(max(q2,1e-300)) if np.isfinite(q2) else np.nan,np.linalg.cond(A)])
    with open(OUT/"covariance_certificate.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["kernel","dimension","n","certificate","assembled_Q_factor","log10_difference","cond2"]); w.writerows(rows)
    print("WROTE",OUT/"covariance_certificate.csv",len(rows),"rows")
if __name__=="__main__": main("--quick" in sys.argv)
