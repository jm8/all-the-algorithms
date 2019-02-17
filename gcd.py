# Binary GCD algorithm
# https://en.wikipedia.org/wiki/Binary_GCD_algorithm

def gcd(u: int, v: int) -> int:
    if u==v: return u
    if u==0: return v
    if v==0: return u
    
    if (~u & 1):
        if v & 1:
            return gcd(u >> 1, v)
        else:
            return gcd(u >> 1, v >> 1) << 1
        
    if ~v & 1:
        return gcd(u, v >> 1)
    
    if (u > v):
        return gcd((u - v) >> 1, v)
    else:
        return gcd((v - u) >> 1, u)
