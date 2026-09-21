import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np, scipy.linalg as la
from ibmi_cert import *

def test_a_only_certificate_matches_sweep():
    rng=np.random.default_rng(7); X=rng.standard_normal((12,12)); A=X.T@X+np.eye(12)
    B=contiguous_overlapping_blocks(12,4,.10)
    g=certificate(A,B); Q=sweep_operator(A,B)
    assert abs(g-max(abs(la.eigvals(Q)))**2) < 1e-10

def test_ibmi_congruence_one_sweep():
    rng=np.random.default_rng(8); X=rng.standard_normal((10,10)); A=X.T@X+2*np.eye(10); H=la.inv(A)
    B=contiguous_overlapping_blocks(10,3,.1); X0=np.eye(10)
    X1=ibmi_sweep(A,X0,B); Q=sweep_operator(A,B)
    assert la.norm((X1-H)-Q@(X0-H)@Q.T) < 1e-10

def test_singleton_cycle_formula():
    H=np.array([[2.,.2,.1],[.2,2.,.3],[.1,.3,2.]])
    A=la.inv(H); B=[np.setdiff1d(np.arange(3),[i]) for i in range(3)]
    assert abs(certificate(A,B)-singleton_cycle_factor(H,(0,1,2))) < 1e-12
