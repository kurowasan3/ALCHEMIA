CODE = 'ORG-ARO-038'
# ---- 共通チェック（g10・難関編） ----
import sys, os, re, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
ITEM = os.path.join('/home/claude/chem/out2/items', CODE)
J = json.load(open(os.path.join(ITEM, 'item.json'), encoding='utf-8'))
FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
assert J['field'] in FIELDS, J['field']
assert J['difficulty'] in (4, 5) and 20 <= J['time_min'] <= 40
fm = J['format'].split(',')
assert 1 <= len(fm) <= 2 and all(f in 'qa equation compare condition pathway calc structure'.split() for f in fm)
assert 3 <= len(J['tags'].split(',')) <= 6
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
for k in ('problem', 'answer', 'solution'):
    assert J[k].count('{') == J[k].count('}'), k
assert J['solution'].rstrip().endswith('\\hfill$\\square$')
print(f'[共通] {CODE}: 反応式 {n_eq} 本つり合いOK／小問 {n_items} 個と答えの番号一致／形式OK')
def formula(smi):
    return CalcMolFormula(Chem.MolFromSmiles(smi))
def mw(smi, table={'H':1.0,'C':12,'N':14,'O':16,'Cl':35.5,'Br':80,'S':32,'Na':23,'K':39}):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    return sum(table[a.GetSymbol()] for a in m.GetAtoms())
def canon(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))
# ---- ここから個別の検証 ----

from rdkit.Chem import AllChem
alk = ['CCCCC','C(C)CCC','CC(C)CC','CCC(C)C','C(CC)CC','C(C)(C)CC','C(C)C(C)C','CC(C)(C)C']
res = set()
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
for a in alk:
    m = Chem.MolFromSmiles('c1ccccc1' + a)
    for b in m.GetBonds():
        if b.GetBeginAtom().GetIsAromatic() or b.GetEndAtom().GetIsAromatic(): continue
        rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.TRIPLE)
        x = rw.GetMol()
        try: Chem.SanitizeMol(x)
        except Exception: continue
        res.add(Chem.MolToSmiles(x))
assert all(formula(s) == 'C11H12' for s in res)
term = Chem.MolFromSmarts('[CH]#C')
internal = sorted(s for s in res if not Chem.MolFromSmiles(s).HasSubstructMatch(term))
assert len(internal) == 5, internal
def branched(s):
    m = Chem.MolFromSmiles(s)
    return any(a.GetDegree() > 2 for a in m.GetAtoms() if not a.GetIsAromatic())
linear = sorted(s for s in internal if not branched(s))
assert len(linear) == 3
# 水の付加：三重結合の各炭素を C=O に
rxn_products = {}
for s in linear:
    m = Chem.MolFromSmiles(s)
    tb = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.TRIPLE][0]
    prods = []
    for ketC, hC in ((tb.GetBeginAtomIdx(), tb.GetEndAtomIdx()), (tb.GetEndAtomIdx(), tb.GetBeginAtomIdx())):
        rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(ketC, hC).SetBondType(Chem.BondType.SINGLE)
        o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(ketC, o, Chem.BondType.DOUBLE)
        x = rw.GetMol(); Chem.SanitizeMol(x)
        aryl_cation = any(n.GetIsAromatic() for n in x.GetAtomWithIdx(ketC).GetNeighbors())
        prods.append((Chem.MolToSmiles(x), aryl_cation))
    rxn_products[s] = prods
mk = Chem.MolFromSmarts('[CH3]C(=O)[#6]')
iodo = {s: any(Chem.MolFromSmiles(p).HasSubstructMatch(mk) for p, _ in v) for s, v in rxn_products.items()}
biased = {s: any(a for _, a in v) for s, v in rxn_products.items()}
ok = [s for s in linear if not iodo[s] and biased[s]]
assert ok == [canon('CCCC#Cc1ccccc1')], ok
assert iodo[canon('CC#CCCc1ccccc1')] and not iodo[canon('CCC#CCc1ccccc1')] and not biased[canon('CCC#CCc1ccccc1')]
P = canon('CCCC#Cc1ccccc1'); X = canon('CCCCC(=O)c1ccccc1'); Y = canon('CCCC(=O)Cc1ccccc1')
pr = dict(rxn_products[P]); assert pr[X] is True and pr[Y] is False
assert formula(X) == formula(Y) == 'C11H14O'
R = Chem.MolFromSmiles('CCCCC(O)c1ccccc1')
assert len(Chem.FindMolChiralCenters(R, includeUnassigned=True)) == 1
S = Chem.MolFromSmiles('CCCC(O)Cc1ccccc1')
assert len(Chem.FindMolChiralCenters(S, includeUnassigned=True)) == 1
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
alkene = Chem.MolFromSmiles('CCCC=Cc1ccccc1')
assert len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(alkene)}) == 2
assert len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(R)}) == 2
for eq in ['C6H5C#CCH2CH2CH3 + H2O -> C6H5COCH2CH2CH2CH3', 'C6H5CH=CHCH2CH2CH3 + H2O -> C6H5CH(OH)CH2CH2CH2CH3']:
    assert check_equation(eq) is None
print('ORG-ARO-038 個別検証OK', internal)
