# THE-BND-095 verification: average valence electrons = 4 rule, EN differences
from fractions import Fraction as F
val={'Al':3,'P':5,'Cd':2,'Te':6,'Ga':3,'S':6,'Cu':1,'In':3,'Zn':2,'Ge':4,'Sn':4,'O':6,'K':1,'Cl':7,'Mg':2,'Si':4,'As':5,'Se':6,'Br':7}
cmp={'AlP':{'Al':1,'P':1},'CdTe':{'Cd':1,'Te':1},'Ga2S3':{'Ga':2,'S':3},'CuInS2':{'Cu':1,'In':1,'S':2},
     'ZnGeP2':{'Zn':1,'Ge':1,'P':2},'SnO2':{'Sn':1,'O':2},'KCl':{'K':1,'Cl':1},'Mg2Si':{'Mg':2,'Si':1}}
avg={k:F(sum(val[e]*n for e,n in v.items()),sum(v.values())) for k,v in cmp.items()}
print(avg)
four=[k for k,v in avg.items() if v==4]; print(four)
assert four==['AlP','CdTe','CuInS2','ZnGeP2','KCl']
EN={'K':0.8,'Cl':3.2,'Ge':2.0,'Ga':1.8,'As':2.2,'Zn':1.7,'Se':2.6,'Cu':1.9,'Br':3.0}
print('KCl diff',EN['Cl']-EN['K'])
# Cu_x Zn Sn S4
import sympy as sp
x=sp.symbols('x'); sol=sp.solve(sp.Eq((x*1+2+4+4*6),4*(x+1+1+4)),x); print('x',sol); assert sol==[2]
d={'Ge':0,'GaAs':round(EN['As']-EN['Ga'],2),'ZnSe':round(EN['Se']-EN['Zn'],2),'CuBr':round(EN['Br']-EN['Cu'],2)}
print(d); order=sorted(d,key=d.get); print(order); assert order==['Ge','GaAs','ZnSe','CuBr']
# literature band gaps (eV, 300 K) for consistency: Ge 0.66, GaAs 1.42, ZnSe 2.7, CuBr 3.0
lit={'Ge':0.66,'GaAs':1.42,'ZnSe':2.7,'CuBr':3.0}; assert sorted(lit,key=lit.get)==order
print('OK')
