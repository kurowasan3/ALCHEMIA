# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/outh/items/THE-RDX-054/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from sympy import symbols, solve, Rational as R
k,m,n,a,f=symbols('k m n a f')
n0=R(2,1000); co2=R(6,1000); h2o=R(6,1000); o2=R(8,1000)
s=solve([n0*k-co2, n0*m-2*h2o, n0*n+2*o2-2*co2-h2o],[k,m,n]); print(s); assert s=={k:3,m:6,n:1}
# 電子数（炭素の酸化数の差）：H=+1,O=-2
def e(C,H,O): return 4*C + H - 2*O   # 完全酸化で失う電子数（C:+4 にする）
assert e(3,6,1)==16 and e(2,4,2)==8 and e(1,2,2)==2
re_=n0-R(8,10000); C=3*re_
sol=solve([2*a+f+R(150,100000)-C, 16*re_-8*a-2*f-4*R(315,100000)],[a,f]); print(sol)
assert sol[a]==R(6,10000) and sol[f]==R(9,10000)
H_w=6*re_-4*sol[a]-2*sol[f]; dw=H_w/2; assert dw==R(15,10000)
O_in=re_+2*R(315,100000); O_out=2*R(150,100000)+2*sol[a]+2*sol[f]+dw; assert O_in==O_out
a1=(C-R(150,100000))/2; o2_1=(16*re_-8*a1)/4; assert a1==R(105,100000) and o2_1==R(270,100000)
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
iso={'アセトン':'CC(C)=O','プロパナール':'CCC=O','アリルアルコール':'C=CCO','1-プロペン-1-オール':'CC=CO','1-プロペン-2-オール':'CC(O)=C',
 'メチルビニルエーテル':'COC=C','シクロプロパノール':'OC1CC1','オキセタン':'C1COC1','メチルオキシラン':'CC1CO1'}
pat1=Chem.MolFromSmarts('[CH3][CX3](=O)'); pat2=Chem.MolFromSmarts('[CH3][CX4H1][OH]')
cans=set(); hits=[]
for nm,smi in iso.items():
    mm=Chem.MolFromSmiles(smi); assert CalcMolFormula(mm)=='C3H6O'; cans.add(Chem.MolToSmiles(mm))
    if mm.HasSubstructMatch(pat1) or mm.HasSubstructMatch(pat2): hits.append(nm)
assert len(cans)==9 and hits==['アセトン'], hits
print('数値・構造 OK')

print('ALL OK')
