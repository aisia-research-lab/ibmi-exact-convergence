import itertools
import numpy as np
import scipy.linalg as la
from .core import cache_local_maps, certificate_from_cache

def singleton_cycle_factor(H,order):
    p=1.0
    for k in range(len(order)):
        p*=abs(H[order[(k+1)%len(order)],order[k]])
    d=np.prod([H[i,i] for i in order])
    return float((p/d)**2)

def two_opt_singleton(H,order):
    order=list(order); cur=singleton_cycle_factor(H,order); changed=True
    while changed:
        changed=False; n=len(order)
        for i in range(n-1):
            for j in range(i+2,n if i else n-1):
                cand=order[:i+1]+list(reversed(order[i+1:j+1]))+order[j+1:]
                val=singleton_cycle_factor(H,cand)
                if val<cur:
                    order,cur=cand,val; changed=True; break
            if changed: break
    return tuple(order),cur

def exhaustive_fixed_block_scores(A,blocks):
    cache=cache_local_maps(A,blocks)
    rows=[]
    for order in itertools.permutations(range(len(blocks))):
        rows.append((order,certificate_from_cache(cache,order)))
    return rows
