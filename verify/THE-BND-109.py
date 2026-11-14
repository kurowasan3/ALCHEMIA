# THE-BND-109 verification: Born-Haber vs Born-Lande lattice energies of NaBr and AgBr
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['Na(s) -> Na(g)','Na(g) -> Na+(g) + e-','Br(g) + e- -> Br-(g)','NaBr(s) -> Na+(g) + Br-(g)',
          'Ag(s) -> Ag(g)','AgBr(s) -> Ag+(g) + Br-(g)','Br2(l) -> 2Br(g)','2Na(s) + Br2(l) -> 2NaBr(s)']:
    r=check_equation(e); print(e,'->',r or 'OK'); assert r is None
sub={'Na':107,'Ag':285}; IE={'Na':496,'Ag':731}; dHf={'Na':-361,'Ag':-100}; Brg=112; EA=325
Uexp={m:-dHf[m]+sub[m]+IE[m]+Brg-EA for m in sub}; print('Uexp',Uexp); assert Uexp=={'Na':751,'Ag':903}
NAke2=1.389e5; M=1.748; a={'Na':597,'Ag':577}; nB={'Na':(7+10)/2,'Ag':(10+10)/2}
EV={m:NAke2*M/(a[m]/2) for m in a}; print('EV',EV)
Ucal={m:EV[m]*(1-1/nB[m]) for m in a}; print('Ucal',Ucal)
assert round(Ucal['Na'])==718 and round(Ucal['Ag'])==757
dev={m:(Uexp[m]-Ucal[m])/Uexp[m]*100 for m in a}; print('dev %',dev)
assert round(dev['Na'],1)==4.4 and round(dev['Ag'])==16
assert EV['Ag']<Uexp['Ag']   # even with no repulsion the ionic model falls short
r_need=NAke2*M*(1-1/10)/Uexp['Ag']; print('r0 needed',r_need); assert round(r_need)==242
print('max (n->inf) shortfall', Uexp['Ag']-EV['Ag'])
print('OK')
