CODE = 'ORG-ARO-028'
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

from rdkit.Chem import AllChem
d = {'A':'c1ccccc1C(O)C(C)C','A2':'c1ccccc1CC(O)(C)C','B':'c1ccccc1C=C(C)C','E':'c1ccccc1CC(C)=C','C':'O=Cc1ccccc1','D':'CC(C)=O'}
f = {k: formula(v) for k, v in d.items()}; print(f)
assert f == {'A':'C10H14O','A2':'C10H14O','B':'C10H12','E':'C10H12','C':'C7H6O','D':'C3H6O'}
cc = lambda s: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
assert cc(d['A']) == 1 and cc(d['A2']) == 0
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
assert len(list(EnumerateStereoisomers(Chem.MolFromSmiles(d['B'])))) == 1
def dehyd(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]; c = o.GetNeighbors()[0]
    for nb in c.GetNeighbors():
        if nb.GetSymbol() == 'C' and not nb.GetIsAromatic() and nb.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
            rw.RemoveAtom(o.GetIdx()); x = rw.GetMol(); Chem.SanitizeMol(x); out.add(Chem.MolToSmiles(x))
    return out
print(dehyd(d['A']), dehyd(d['A2']))
assert dehyd(d['A']) == {canon(d['B'])} and dehyd(d['A2']) == {canon(d['B']), canon(d['E'])}
rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]=O.O=[C:2]')
def ozo(s):
    ps = rxn.RunReactants((Chem.MolFromSmiles(s),))
    out = set()
    for p in ps:
        for x in p:
            Chem.SanitizeMol(x)
        out.add(tuple(sorted(Chem.MolToSmiles(x) for x in p)))
    return out
print(ozo(d['B']), ozo(d['E']))
assert ozo(d['B']) == {tuple(sorted([canon('CC(C)=O'), canon('O=Cc1ccccc1')]))}
assert ozo(d['E']) == {tuple(sorted([canon('C=O'), canon('CC(=O)Cc1ccccc1')]))}
print('ORG-ARO-028 個別検証OK')
