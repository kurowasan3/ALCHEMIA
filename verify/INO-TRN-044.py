# INO-TRN-044 鉄の酸化物の混合物：組成は決まらないが Fe(II)/Fe(III) の量は決まる
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, symbols, solve
Fe,O=56,16
m=R(600,100); mox=R(640,100)
nFe=mox/(2*Fe+3*O)*2
mFe=nFe*Fe
nO=(m-mFe)/O
print('Fe mol',nFe,'Fe g',float(mFe),'O mol',float(nO))
n3=2*(nO-nFe); n2=nFe-n3
print('Fe3+',float(n3),'Fe2+',float(n2), 'charge', 2*n2+3*n3, 2*nO)
# 組成の自由度：FeO a, Fe3O4 b, Fe2O3 c
a,b,c=symbols('a b c')
sol=solve([a+3*b+2*c-nFe, a+4*b+3*c-nO],[a,c],dict=True)[0]
print(sol)
# b の範囲：a>=0, c>=0
import sympy
bmax=min(sympy.solve(sol[a],b)[0], sympy.solve(sol[c],b)[0]) if True else None
print('a=0 at b=',sympy.solve(sol[a],b),' c=0 at b=',sympy.solve(sol[c],b))
for bb in [0,R(1,100),R(15,1000)]:
    aa=sol[a].subs(b,bb); cc=sol[c].subs(b,bb)
    print('b',bb,'a',aa,'c',cc,'mass',72*aa+232*bb+160*cc)
for e in ['FeO + H2 -> Fe + H2O','Fe3O4 + 4H2 -> 3Fe + 4H2O','Fe2O3 + 3H2 -> 2Fe + 3H2O','4FeO + O2 -> 2Fe2O3','4Fe3O4 + O2 -> 6Fe2O3']:
    print(e, check_equation(e))
print('H2 L at STP', float(nO*R(224,10)))
