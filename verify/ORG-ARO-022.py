CODE = 'ORG-ARO-022'
# ---- 共通チェック（g10） ----
import sys, os, re, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
ITEM = os.path.join('/home/claude/chem/out/items', CODE)
J = json.load(open(os.path.join(ITEM, 'item.json'), encoding='utf-8'))
FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
assert J['field'] in FIELDS, J['field']
assert 1 <= J['difficulty'] <= 3
assert all(f in 'qa equation compare condition pathway calc structure'.split() for f in J['format'].split(','))
def _ce_all(s):
    out = []; i = 0
    while True:
        k = s.find('\\ce{', i)
        if k < 0: return out
        d = 0; j = k + 3
        while True:
            if s[j] == '{': d += 1
            elif s[j] == '}':
                d -= 1
                if d == 0: break
            j += 1
        out.append(s[k+4:j]); i = j
n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for e in _ce_all(J[key]):
        arrows = re.findall(r'<=>|->|<-', e)
        if len(arrows) != 1 or '\\boxed' in e or re.search(r'\bR\d?\b|\bn\b|\bAr\b|\bX\b', e) or '...' in e:
            continue
        r = check_equation(e)
        assert r is None, (key, e, r)
        n_eq += 1
n_items = J['problem'].count('\\item')
nums = re.findall(r'(?:^|\s)\((\d+)\)', J['answer'])
if n_items:
    assert [int(x) for x in nums] == list(range(1, n_items + 1)), (n_items, nums)
    assert J['problem'].count('\\begin{enumerate}') == 1
else:
    assert not nums, nums
for k in ('problem_fig', 'solution_fig'):
    if J[k]:
        sv = open(os.path.join(ITEM, J[k]), encoding='utf-8').read()
        assert 'viewBox' in sv and '<script' not in sv and 'foreignObject' not in sv and '<image' not in sv
        assert not re.search(r'\son[a-z]+=', sv) and 'href=' not in sv
print(f'[共通] {CODE}: 反応式 {n_eq} 本つり合いOK／小問 {n_items} 個と答えの番号一致／形式OK')
def formula(smi):
    m = Chem.MolFromSmiles(smi); return CalcMolFormula(m)
def mw(smi, table={'H':1.0,'C':12,'N':14,'O':16,'Cl':35.5,'Br':80,'S':32,'Na':23,'K':39}):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    return sum(table[a.GetSymbol()] for a in m.GetAtoms())
def canon(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))
# ---- ここから個別の検証 ----

from sympy import Rational as Q
O = Q(10000-7962-616-664,100)
r = [Q(7962,100)/12, Q(616,100)/1, Q(664,100)/14, O/16]
base = r[2]
print([float(x/base) for x in r])
assert [round(float(x/base)) for x in r] == [14,13,1,1]
A = 'Cc1ccc(NC(=O)c2ccccc2)cc1'
assert formula(A) == 'C14H13NO' and mw(A) == 211
assert formula('OC(=O)c1ccccc1') == 'C7H6O2' and formula('Cc1ccc(N)cc1') == 'C7H9N'
def ring_cl(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for a in m.GetAtoms():
        if a.GetIsAromatic() and a.GetTotalNumHs() == 1:
            rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(17)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
            x = rw.GetMol(); Chem.SanitizeMol(x); out.add(Chem.MolToSmiles(x))
    return len(out)
res = [ring_cl(s) for s in ('Cc1ccccc1N', 'Cc1cccc(N)c1', 'Cc1ccc(N)cc1')]
print(res); assert res == [4,4,2]
print('ORG-ARO-022 個別検証OK')
