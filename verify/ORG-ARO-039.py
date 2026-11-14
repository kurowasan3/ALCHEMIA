CODE = 'ORG-ARO-039'
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

eqs = ['C6H5C(CH3)2OOH -> C6H5COCH3 + CH3OH', 'C6H5C(CH3)2OOH -> C6H5OH + CH3COCH3',
       'C6H5C(CH3)(C2H5)OOH -> C6H5OH + CH3COC2H5', 'CH3C6H4C(CH3)2OOH -> CH3C6H4OH + CH3COCH3',
       'C6H4(C(CH3)2OOH)2 -> C6H4(OH)2 + 2CH3COCH3']
for e in eqs: assert check_equation(e) is None, e
benzylic_one_H = Chem.MolFromSmarts('c[CX4H1]')
mono = {'butyl': 'CCCCc1ccccc1', 'isobutyl': 'CC(C)Cc1ccccc1', 'sec-butyl': 'CCC(C)c1ccccc1', 'tert-butyl': 'CC(C)(C)c1ccccc1'}
assert all(formula(s) == 'C10H14' for s in mono.values()) and len({canon(s) for s in mono.values()}) == 4
hit = [k for k, s in mono.items() if Chem.MolFromSmiles(s).HasSubstructMatch(benzylic_one_H)]
assert hit == ['sec-butyl'], hit
para = {'p-diethyl': 'CCc1ccc(CC)cc1', 'p-propyltoluene': 'CCCc1ccc(C)cc1', 'p-cymene': 'Cc1ccc(C(C)C)cc1'}
assert all(formula(s) == 'C10H14' for s in para.values())
hit = [k for k, s in para.items() if Chem.MolFromSmiles(s).HasSubstructMatch(benzylic_one_H)]
assert hit == ['p-cymene']
assert formula('Oc1ccc(C)cc1') == 'C7H8O'
assert mw('CC(C)c1cccc(C(C)C)c1') == 162 and mw('Oc1cccc(O)c1') == 110
m = 0.100 * 0.95**2 * 110; assert round(m, 1) == 9.9, m
assert abs(0.95**2 + 2*0.95*0.05 + 0.05**2 - 1) < 1e-12
assert formula('CC(=O)c1cccc(O)c1') == 'C8H8O2' and formula('CC(=O)c1cccc(C(C)=O)c1') == 'C10H10O2'
print('ORG-ARO-039 個別検証OK', m)
