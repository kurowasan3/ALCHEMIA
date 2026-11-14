CODE = 'ORG-ARO-037'
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

from fractions import Fraction
from itertools import combinations
for eq in ['MnO4- + 2H2O + 3e- -> MnO2 + 4OH-', 'MnO4- + 8H+ + 5e- -> Mn^2+ + 4H2O']:
    assert check_equation(eq) is None
# 側鎖酸化の電子数：側鎖炭素の酸化数和（H は +1, C-C は 0, O は -2）
def side_ox(m, side_atoms):
    s = 0
    for i in side_atoms:
        a = m.GetAtomWithIdx(i)
        s -= a.GetTotalNumHs()
        for nb in a.GetNeighbors():
            if nb.GetSymbol() == 'O':
                b = m.GetBondBetweenAtoms(i, nb.GetIdx()); s += int(b.GetBondTypeAsDouble())
    return s
alkyls = {1:['C'],2:['CC'],3:['CCC','C(C)C'],4:['CCCC','C(C)CC','CC(C)C','C(C)(C)C'],5:['CCCCC','C(C)CCC','CC(C)CC','CCC(C)C','C(CC)CC','C(C)(C)CC','C(C)C(C)C','CC(C)(C)C']}
for n, lst in alkyls.items():
    for a in lst:
        m = Chem.MolFromSmiles('c1ccccc1' + a)
        side = [x.GetIdx() for x in m.GetAtoms() if not x.GetIsAromatic()]
        benz = m.GetAtomWithIdx(6)
        if benz.GetTotalNumHs() == 0: continue
        before = side_ox(m, side)
        after = 3 + 4 * (n - 1)
        assert after - before == 6 * n, (a, after - before)
print('電子数 6n OK')
sol = [(k, Fraction(67 * k - 78, 14)) for k in range(1, 7)]
ok = [(k, N) for k, N in sol if N >= 0 and N.denominator == 1]
assert ok == [(2, 4)], ok
# 質量の整合
e = Fraction(400, 10000) * 3
assert e == Fraction(12, 100) and e / 6 == Fraction(2, 100)
assert abs(1.62 / 0.0100 - 162) < 1e-9 and abs(2.22 / 0.0100 - 222) < 1e-9
# 位置異性体
iso = set()
import itertools
for pos in itertools.permutations(range(6), 3):
    sub = ['', '', '', '', '', '']
    sub[pos[0]] = '(C)'; sub[pos[1]] = '(C)'; sub[pos[2]] = '(C(C)(C)C)'
    smi = 'c1' + ''.join(f'c{sub[i]}' for i in range(1, 6))
    smi = 'c' + sub[0] + '1' + ''.join(f'c{sub[i]}' for i in range(1, 6)) + '1'
    iso.add(canon(smi))
assert len(iso) == 6
assert all(formula(s) == 'C12H18' and mw(s) == 162 for s in iso)
patt = Chem.MolFromSmarts('[CH3]c:c[CH3]')
rest = [s for s in iso if not Chem.MolFromSmiles(s).HasSubstructMatch(patt)]
assert len(rest) == 4, rest
acid = Chem.MolFromSmiles('OC(=O)c1cc(C(=O)O)cc(C(C)(C)C)c1')
assert mw('OC(=O)c1cc(C(=O)O)cc(C(C)(C)C)c1') == 222
assert abs(0.120 / 5 - 0.0240) < 1e-12
print('ORG-ARO-037 個別検証OK', sorted(rest))
