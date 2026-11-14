"""THE-ELC-041 検証：電気分解の結果から水溶液を決める・極性の反転（SymPy／chemlib）"""
import sys
from sympy import Rational as R, log, N
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

ne = R(4, 100)            # 各容器に流した電子 0.0400 mol
c0 = R(1, 10)             # 各水溶液 0.100 mol/L, 1.00 L
Ag, Cu, Vm = 108, R(635, 10), R(224, 10)
# 候補ごとの生成物（陰極, 陽極）と電解後の液性
cand = {
 '硝酸銀':   dict(cath='Ag', an='O2', H=+ne),          # 4Ag+ + 2H2O -> 4Ag + O2 + 4H+
 '塩化銅(II)': dict(cath='Cu', an='Cl2', H=0),
 '硫酸ナトリウム': dict(cath='H2', an='O2', H=0),
 '塩酸':     dict(cath='H2', an='Cl2', H=-ne),
 '塩化ナトリウム': dict(cath='H2', an='Cl2', H=-ne),   # OH- 生成（塩基性）
}
initH = {'硝酸銀':0,'塩化銅(II)':0,'硫酸ナトリウム':0,'塩酸':c0,'塩化ナトリウム':0}
# 条件 (a) A,C 陰極の質量増加 A>C　(b) B,D 陰極で同じ気体　(c) C,D 陽極 Cl2，A,B 陽極 O2　(d) 混合後 B 中性・D 酸性
def mass(s): return {'Ag': ne*Ag, 'Cu': ne/2*Cu}.get(cand[s]['cath'], 0)
sols = []
names = list(cand)
import itertools
for A,B,C,D in itertools.permutations(names, 4):
    if not (mass(A)>0 and mass(C)>0 and mass(A)>mass(C)): continue
    if not (cand[B]['cath']=='H2' and cand[D]['cath']=='H2'): continue
    if not (cand[C]['an']=='Cl2' and cand[D]['an']=='Cl2' and cand[A]['an']=='O2' and cand[B]['an']=='O2'): continue
    hB = initH[B]+cand[B]['H']; hD = initH[D]+cand[D]['H']
    if not (hB==0 and hD>0): continue
    sols.append((A,B,C,D))
print('解', sols); assert sols == [('硝酸銀','硫酸ナトリウム','塩化銅(II)','塩酸')]
mA, mC = ne*Ag, ne/2*Cu
print('A 陰極', float(mA), 'C 陰極', float(mC)); assert round(float(mA),2)==4.32 and round(float(mC),2)==1.27
VD = (ne/2 + ne/2)*Vm; print('D 気体合計 L', float(VD)); assert VD == R(896,1000)
pHA = -log(ne, 10); print('A pH', float(N(pHA))); assert round(float(N(pHA)),2) == 1.40
pHD = -log(c0-ne, 10); print('D pH (参考)', float(N(pHD)))
# 極性反転：A で Ag が全部溶けるまでの電子 = 析出した Ag の物質量
print('反転後 Ag が溶け切るまでの電子', ne)
for eq in ['Ag+ + e- -> Ag','2H2O -> O2 + 4H+ + 4e-','Cu^2+ + 2e- -> Cu','2Cl- -> Cl2 + 2e-',
           '2H+ + 2e- -> H2','2H2O + 2e- -> H2 + 2OH-','Ag -> Ag+ + e-','Cu -> Cu^2+ + 2e-',
           '4AgNO3 + 2H2O -> 4Ag + O2 + 4HNO3']:
    r = check_equation(eq); print(eq, r or 'OK'); assert r is None
print('THE-ELC-041 all OK')
