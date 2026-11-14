CODE = 'THE-EQL-049'
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

import sympy as sp
for e in ['C6H6 + C2H4 -> C6H5C2H5', 'C6H5C2H5 + C2H4 -> C6H4(C2H5)2']:
    assert check_equation(e) is None
x = sp.symbols('x', positive=True)          # x = n_B / n0
M = sp.Function('M')
sol = sp.dsolve(sp.Eq(M(x).diff(x), -1 + 2*M(x)/x), M(x), ics={M(1): 0})
assert sp.simplify(sol.rhs - x*(1 - x)) == 0
f = sp.symbols('f')
D = 1 - f - f*(1 - f)
assert sp.expand(D - (1 - f)**2) == 0
Euse = sp.expand(f*(1 - f) + 2*(1 - f)**2)
assert sp.factor(Euse - (1 - f)*(2 - f)) == 0
fs = [s for s in sp.solve(sp.Eq(Euse, sp.Rational(1, 10)), f) if 0 < s <= 1][0]
ratio = float((1 - fs)/fs); assert round(ratio, 2) == 0.10, ratio
share = float(2*(1 - fs)**2/sp.Rational(1, 10)); assert abs(share - 0.168) < 0.002
f4 = sp.Rational(100, 101); E4 = float((1 - f4)*(2 - f4)); assert round(E4, 3) == 0.010 and E4 <= 0.010
# 数値積分による独立確認（E を少しずつ加える）
def sim(Etot, steps=200000):
    B, Mm, Dd = 1.0, 0.0, 0.0
    dE = Etot/steps
    for _ in range(steps):
        a = B/(B + 2*Mm); B -= dE*a; Mm += dE*a - dE*(1 - a); Dd += dE*(1 - a)
    return B, Mm, Dd
B, Mm, Dd = sim(0.1)
assert abs(Dd/Mm - ratio) < 1e-3 and abs(B - float(fs)) < 1e-4
print('THE-EQL-049 個別検証OK', float(fs), ratio, share, E4, Dd/Mm)
