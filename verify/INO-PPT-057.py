# INO-PPT-057 verification: high/low spin by brute-force minimisation
def configs(n):
    for x in range(0,7):
        y=n-x
        if 0<=y<=4: yield x,y
def pairs(k,orbs): return max(0,k-orbs)
def best(n,D,P):
    res=[]
    for x,y in configs(n):
        E=-0.4*D*x+0.6*D*y+(pairs(x,3)+pairs(y,2))*P
        unp=(min(x,3)-pairs(x,3))+(min(y,2)-pairs(y,2))
        res.append((round(E,6),unp,x,y))
    return sorted(res)
for n in (4,5,6,7):
    hs=[r for r in best(n,1,0) if True]
    # HS = max unpaired config, LS = max x
    cs=best(n,100,50)
    hsE=max(cs,key=lambda r:(r[1],-r[0])); lsE=max(cs,key=lambda r:r[2])
    D,P=100,50
    print(n,'LS-HS',lsE[0]-hsE[0],'(P-D)=',P-D)
tab={'Cr(NH3)6':(3,258,282,3),'Mn(H2O)6':(5,100,304,5),'Fe(H2O)6':(6,124,210,4),'Fe(CN)6':(6,395,210,0),'CoF6':(6,155,251,4),'Co(NH3)6':(6,274,251,0)}
for k,(n,D,P,exp) in tab.items():
    b=best(n,D,P)[0]; print(k,b); assert b[1]==exp, k
for n in (3,8):
    print(n,'distinct unpaired among minima over Delta/P', {best(n,D,P)[0][1] for D in (50,500) for P in (50,500)})
lam=6.63e-34*3.00e8*6.02e23/210e3*1e9; print('lambda0',lam); assert round(lam)==570
# Ni(II) d8 tetrahedral: low 2 orbitals, high 3, gap << P
def tet(n,D,P):
    r=[]
    for x in range(0,5):
        y=n-x
        if 0<=y<=6:
            E=D*y+(max(0,x-2)+max(0,y-3))*P; u=(min(x,2)-max(0,x-2))+(min(y,3)-max(0,y-3)); r.append((E,u))
    return min(r)
def sq(n,D,P):
    r=[]
    for x in range(0,9):
        y=n-x
        if 0<=y<=2:
            E=D*y+(max(0,x-4)+max(0,y-1))*P; u=(min(x,4)-max(0,x-4))+(min(y,1)-max(0,y-1)); r.append((E,u))
    return min(r)
print('tet',tet(8,10,200),'sq',sq(8,1000,200)); assert tet(8,10,200)[1]==2 and sq(8,1000,200)[1]==0
print('OK')
