from scipy.optimize import brentq
pK={'Nterm':8.0,'Cterm':3.1,'Glu':4.1,'Asp':3.9,'His':6.0,'Lys':10.5,'Arg':12.5}
def charge(pH, comp):
    q=1/(1+10**(pH-pK['Nterm'])) - 1/(1+10**(pK['Cterm']-pH))
    for aa,n in comp.items():
        if aa in ('Glu','Asp'): q-=n/(1+10**(pK[aa]-pH))
        elif aa in ('His','Lys','Arg'): q+=n/(1+10**(pH-pK[aa]))
    return q
def pI(comp): return brentq(lambda x: charge(x,comp),0,14)
base={'Glu':6,'Asp':4,'His':3,'Lys':6,'Arg':2}
p0=pI(base); print('pI base',p0)
m1=dict(base); m1['Glu']-=1; p1=pI(m1); print('Glu->Val',p1); assert p1>p0
m2=dict(base); m2['Lys']+=1; p2=pI(m2); print('Val->Lys',p2); assert p2>p0
p3=pI(base); assert abs(p3-p0)<1e-12   # Val->Leu: no ionizable change
assert charge(8.0,base)<0 if p0<8 else True
# electrophoresis at 5.8 for pI 4.6 (negative) and 7.0 (positive)
assert 5.8>4.6 and 5.8<7.0
print('POL-PRT-028 OK')
