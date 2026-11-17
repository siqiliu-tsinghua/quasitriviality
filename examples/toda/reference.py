"""Exact CP1 Toda and DZ reference identities in translation-invariant Riemann jets."""
import sympy as sp

R = sp.Rational
d=sp.Symbol("d", positive=True)  # r-s=4exp(u/2)>0 on the real chart
N=14
p,q=(sp.symbols(f"{prefix}1:{N+1}") for prefix in ("p","q"))
p1,p2,p3,p4=p[:4]
q1,q2,q3,q4=q[:4]
h=d/4
S=sp.S

def can(f):
    """Here all denominators should be Laurent monomials in d,p1,q1."""
    return sp.expand(f)

def jet_action(f,dot_d,dot_p,dot_q):
    fs=f.free_symbols
    if fs.intersection({p[-1],q[-1]}):
        raise ValueError("reference verifier needs more jet symbols")
    result=dot_d*f.diff(d) if d in fs else S.Zero
    for var,vel in zip(p+q, dot_p+dot_q):
        if var in fs:
            result+=vel*f.diff(var)
    return can(result)

def D(f,n=1):
    for _ in range(n):
        f=jet_action(f,p1-q1,p[1:]+(0,),q[1:]+(0,))
    return f

rv=[D(h*p1)]
sv=[D(-h*q1)]
for k in range(1,N-1):
    rv.append(D(rv[-1]));sv.append(D(sv[-1]))
rv.append(0)
sv.append(0)
def E(f,n=1):
    for _ in range(n):
        f=jet_action(f,h*(p1+q1),rv,sv)
    return f

# DZ 1998, (1.16) and Example 4 G-function.
F1=(sp.log(p1)+sp.log(q1)-2*sp.log(d/4))/24
# DZ 2001, Eq. (3.10.116); symmetric under (p,q,d) -> (q,p,-d).
F2=(4*p2**3*d/(5*p1**4)-4*q2**3*d/(5*q1**4)-p2*q2/(4*p1*q1)
    +3*p2*(p2*q1/2-R(7,5)*p3*d)/(4*p1**3)
    +3*q2*(q2*p1/2+R(7,5)*q3*d)/(4*q1**3)
    +(R(33,10)*p2**2-R(9,10)*p3*q1+p2*q2/10+p4*d)/(4*p1**2)
    +(R(33,10)*q2**2-R(9,10)*q3*p1+p2*q2/10-q4*d)/(4*q1**2)
    -(R(17,5)*p3+q3/2)/(4*p1)-(R(17,5)*q3+p3/2)/(4*q1)
    -(p1**3/q1+q1**3/p1)/(10*d**2)
    -(p1**2-R(11,5)*p1*q1+q1**2)/(d**2)
    +(p2-q2)*(q1/(5*p1)+p1/(5*q1)+1)/d)/24**2
F2=can(F2)

# A=sinh(εD/2)/(εD/2) of the centered Toda map (formal operator).
A=(S.One,R(1,24),R(1,1920),R(1,322560))
A2=(S.One,R(1,12),R(1,360),R(1,20160))
B=(S.One,R(-1,24),R(7,5760),R(-31,967680))
B2=(S.One,R(-1,12),R(1,240),R(-1,6048))
v=[None,D(E(F1)),D(E(F2))]
u=[None,D(F1,2),D(F2,2)]
vder = lambda j: (p[j-1]+q[j-1])/2
uder = lambda j: D(2*(p1-q1)/d,j-1)
U2=can(u[1]+A2[1]*uder(2))
U4=can(u[2]+A2[1]*D(u[1],2)+A2[2]*uder(4))
V2=can(v[1]+A[1]*vder(2))
V4=can(v[2]+A[1]*D(v[1],2)+A[2]*vder(4))

def reference_residuals():
    # Exact tau identity E² F_g = [ε^{2g}] exp(U), with the leading
    # E² F_0 = e^u=(d/4)². No integration needed for this verification.
    return can(E(F1,2)-h**2*U2),can(E(F2,2)-h**2*(U4+U2**2/2))

def canon_from_centered(U,V):
    # Inverse quasi-Miura, i.e. dispersionless -> centered Toda.
    return (can(V[0]+h*U[0]),can(V[0]-h*U[0])), (can(V[1]+h*(U[1]+U[0]**2/4)),can(V[1]-h*(U[1]+U[0]**2/4)))

if __name__ == '__main__':
    print("terms F2, U2, U4",*[len(sp.Add.make_args(a)) for a in (F2,U2,U4)],flush=True)
    print("reference Toda tau residuals eps2/eps4:",reference_residuals(),flush=True)
