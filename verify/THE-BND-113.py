# THE-BND-113 verification: structure factor of rock-salt MX, isoelectronic extinction, KCl identification
import itertools, math
X = [(0,0,0),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5)]
M = [(0.5,0,0),(0,0.5,0),(0,0,0.5),(0.5,0.5,0.5)]
def F(h,k,l,fM,fX):
    c = lambda p: math.cos(2*math.pi*(h*p[0]+k*p[1]+l*p[2]))
    return fX*sum(c(p) for p in X) + fM*sum(c(p) for p in M)
fM, fX = 7.0, 13.0
for h,k,l in itertools.product(range(7), repeat=3):
    v = F(h,k,l,fM,fX); par = {h%2,k%2,l%2}
    exp = 4*(fM+fX) if par == {0} else (4*(fX-fM) if par == {1} else 0.0)
    assert abs(v-exp) < 1e-9, (h,k,l,v,exp)
# (2) NaCl
r = (F(1,1,1,10,18)/F(2,0,0,10,18))**2; print('NaCl ratio', r); assert round(r,3) == 0.082
# (3) isoelectronic
Z = {'Na':11,'K':19,'Rb':37,'F':9,'Cl':17,'Br':35}
pairs = [('Na','F'),('Na','Cl'),('K','F'),('K','Cl'),('K','Br'),('Rb','Br')]
iso = [m+x for m,x in pairs if Z[m]-1 == Z[x]+1]; print('iso', iso)
assert iso == ['NaF','KCl','RbBr']
# (4) pattern
allowed = sorted({h*h+k*k+l*l for h,k,l in itertools.product(range(8),repeat=3)
                  if (h,k,l)!=(0,0,0) and len({h%2,k%2,l%2})==1})
even = sorted({h*h+k*k+l*l for h,k,l in itertools.product(range(0,9,2),repeat=3) if (h,k,l)!=(0,0,0)})
print('allowed', allowed[:9], 'even', even[:6])
obs = [0.0599,0.1199,0.1798,0.2398,0.2997,0.3597]
rat = [o/obs[0] for o in obs]; print([round(x,3) for x in rat])
assert all(abs(x-(i+1))<0.01 for i,x in enumerate(rat))
assert [s/allowed[0] for s in allowed[:3]] != [1,2,3]          # with odd lines: 3:4:8
assert even[:6] == [4,8,12,16,20,24]
lam = 154.0
a = lam*2/(2*math.sqrt(obs[0])); print('a', a); assert round(a) == 629
Mw = 1.98*6.02e23*(a*1e-10)**3/4; print('M', Mw)
mass = {'NaF':23.0+19.0,'KCl':39.1+35.5,'RbBr':85.5+79.9}
best = min(mass, key=lambda k: abs(mass[k]-Mw)); print(best, mass)
assert best == 'KCl' and abs(Mw-74.6) < 1.0
ap = a/2; m1 = 1.98*6.02e23*(ap*1e-10)**3; print("a'", ap, 'mass per point', m1)
assert round(m1,1) == 37.1
assert all(28 != x*x+y*y+z*z for x,y,z in itertools.product(range(6),repeat=3))
print('OK')
