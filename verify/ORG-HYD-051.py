# 検証スクリプト（x1 g07）  python3 -I ORG-HYD-051.py
ITEM = '/home/claude/chemx/x1/out/items/ORG-HYD-051/item.json'
import sys, json, re
from math import log10, floor, log, exp, sqrt
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
assert it['field'] in FIELDS, it['field']
assert set(it['format'].split(',')) <= {'qa','equation','compare','condition','pathway','calc','structure'}

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

neq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            r = check_equation(b)
            assert r is None, (b, r)
            neq += 1
print('反応式チェック:', neq, '本 OK')
nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans, (nprob, nans)
print('小問', nprob, '答え', nans)

def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def near(a, b, rel=0.01):
    return abs(a - b) <= rel * abs(b)


from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors as D
cand={'P':'C#CCC','Q':'CC#CC','R':'C1=CCC1','S':'CC1=CC1','T':'CC1C=C1','U':'C=C1CC1'}
can=set()
for k,s in cand.items():
    m=Chem.MolFromSmiles(s); assert D.CalcMolFormula(m)=='C4H6',k; can.add(Chem.MolToSmiles(m))
assert len(can)==6
for s in ['C=C=CC','C=CC=C','C1C2CC12']:
    assert D.CalcMolFormula(Chem.MolFromSmiles(s))=='C4H6'
def ozon(s):
    m=Chem.MolFromSmiles(s); out=set()
    for b in m.GetBonds():
        if b.GetBondType()==Chem.BondType.DOUBLE and b.GetBeginAtom().GetSymbol()=='C' and b.GetEndAtom().GetSymbol()=='C':
            rw=Chem.RWMol(m); a1,a2=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
            rw.RemoveBond(a1,a2)
            for a in (a1,a2):
                o=rw.AddAtom(Chem.Atom(8)); rw.AddBond(a,o,Chem.BondType.DOUBLE)
            mm=rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return out
can_=lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
assert ozon(cand['R'])=={can_('O=CCCC=O')}, ozon(cand['R'])
assert ozon(cand['S'])=={can_('CC(=O)CC=O')}
assert ozon(cand['T'])=={can_('CC(C=O)C=O')}
assert ozon(cand['U'])=={can_('O=C1CC1.C=O')}, ozon(cand['U'])
pat=Chem.MolFromSmarts('[CH3]C(=O)')
assert Chem.MolFromSmiles('CC(=O)CC=O').HasSubstructMatch(pat) and not Chem.MolFromSmiles('CC(C=O)C=O').HasSubstructMatch(pat)
hyd=lambda s: can_(s.replace('=','').replace('#',''))
assert hyd(cand['S'])==hyd(cand['T'])==hyd(cand['U'])==can_('CC1CC1') and hyd(cand['R'])!=hyd(cand['S'])

print('ALL OK')
