import numpy as np
import scipy.linalg as la

def complement(n, I):
    return np.setdiff1d(np.arange(n), np.asarray(I,dtype=int), assume_unique=True)

def local_graph_map(A, I):
    """Theorem 5.1: T_I = E_I(-A_II^{-1}A_IJ)+E_J, without forming A^{-1}."""
    A=np.asarray(A,float); n=A.shape[0]; I=np.asarray(I,dtype=int); J=complement(n,I)
    X=la.solve(A[np.ix_(I,I)],A[np.ix_(I,J)],assume_a="pos",check_finite=False)
    T=np.zeros((n,len(J))); T[I,:]=-X; T[J,:]=np.eye(len(J))
    return J,T

def cache_local_maps(A, blocks):
    return [local_graph_map(A,I) for I in blocks]

def return_operator_from_cache(cache, order):
    """Theorem 5.2: exact A-only cyclic return operator M_A."""
    order=tuple(order)
    J0,T0=cache[order[0]]
    P=T0[cache[order[1]][0],:]
    for z in range(1,len(order)-1):
        _,T=cache[order[z]]
        P=T[cache[order[z+1]][0],:]@P
    _,Tlast=cache[order[-1]]
    return Tlast[J0,:]@P

def certificate_from_cache(cache, order):
    M=return_operator_from_cache(cache,order)
    return float(np.max(np.abs(la.eigvals(M)))**2)

def certificate(A, blocks, order=None):
    if order is None: order=tuple(range(len(blocks)))
    return certificate_from_cache(cache_local_maps(A,blocks),order)

def sweep_operator(A, blocks, order=None):
    """Independently assemble Q from local graph maps; used only for validation."""
    if order is None: order=tuple(range(len(blocks)))
    n=A.shape[0]; Q=np.eye(n)
    for k in order:
        J,T=local_graph_map(A,blocks[k])
        Qi=np.zeros((n,n)); Qi[:,J]=T
        Q=Qi@Q
    return Q

def reconstruct_block(A,I,XJ):
    n=A.shape[0]; I=np.asarray(I,dtype=int); J=complement(n,I)
    AI=A[np.ix_(I,I)]
    R=la.solve(AI,A[np.ix_(I,J)],assume_a="pos",check_finite=False)
    X=np.zeros((n,n))
    X[np.ix_(I,I)]=la.inv(AI)+R@XJ@R.T
    X[np.ix_(I,J)]=-R@XJ
    X[np.ix_(J,I)]=X[np.ix_(I,J)].T
    X[np.ix_(J,J)]=XJ
    return (X+X.T)/2

def ibmi_sweep(A,X,blocks,order=None):
    if order is None: order=tuple(range(len(blocks)))
    for k in order:
        I=blocks[k]; J=complement(A.shape[0],I)
        X=reconstruct_block(A,I,X[np.ix_(J,J)])
    return X

def contiguous_overlapping_blocks(n,K=4,overlap=0.05):
    cuts=np.linspace(0,n,K+1,dtype=int); pad=int(round(overlap*n/2))
    return [np.arange(max(0,cuts[k]-pad),min(n,cuts[k+1]+pad)) for k in range(K)]
