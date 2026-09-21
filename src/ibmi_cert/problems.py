import numpy as np
import scipy.linalg as la

def points_1d(n): return np.linspace(0,n**0.9,n)[:,None]
def points_2d(n):
    m=int(round(np.sqrt(n)))
    if m*m!=n: raise ValueError("2D grid requires n to be a perfect square")
    x=np.linspace(0,1,m); X,Y=np.meshgrid(x,x,indexing="ij")
    return np.c_[X.ravel(),Y.ravel()]

def covariance_matrix(points,kind,sigma=.5,tau=3.):
    d=la.norm(points[:,None,:]-points[None,:,:],axis=2)
    kind=kind.upper()
    if kind=="EXP": A=np.exp(-d/5)
    elif kind=="RBF": A=np.exp(-(d*d)/(2*sigma*sigma))
    elif kind=="IQUAD": A=1/np.sqrt(1+d*d)
    elif kind=="M32":
        z=np.sqrt(3)*d/tau; A=(1+z)*np.exp(-z)
    elif kind=="M52":
        z=np.sqrt(5)*d/tau; A=(1+z+z*z/3)*np.exp(-z)
    else: raise ValueError(kind)
    return (A+A.T)/2

def scale_diagonal_congruence(A):
    A=np.asarray(A,float); d=np.sqrt(np.diag(A))
    if np.any(d<=0): raise ValueError("positive diagonal required")
    return A/d[:,None]/d[None,:]
