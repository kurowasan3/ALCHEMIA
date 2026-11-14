CODE = 'INO-AMP-014'
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

import sympy as sp, math
for eq in ['C6H5NO2 + 7H+ + 6e- -> C6H5NH3+ + 2H2O','2H+ + 2e- -> H2','Sn^2+ + 2OH- -> Sn(OH)2','Sn^4+ + 4OH- -> Sn(OH)4',
           'C6H5NH3+ + OH- -> C6H5NH2 + H2O','Sn(OH)2 + 2OH- -> [Sn(OH)4]^2-','Sn(OH)4 + 2OH- -> [Sn(OH)6]^2-',
           '2C6H5NO2 + 3Sn + 14HCl -> 2C6H5NH3Cl + 3SnCl4 + 4H2O']:
    assert check_equation(eq) is None, eq
a, b = sp.symbols('a b')
s = sp.solve([a + b - sp.Rational(4, 100), 2*a + 4*b - (sp.Rational(2, 100)*6 + sp.Rational(5, 1000)*2)], [a, b])
assert s[a] == sp.Rational(15, 1000) and s[b] == sp.Rational(25, 1000)
def stages(a, b, HCl=sp.Rational(3, 10), n=sp.Rational(2, 100)):
    H = HCl - (2*a + 4*b + n)
    s1 = H; s2 = s1 + 2*a + 4*b; s3 = s2 + n; s4 = s3 + 2*(a + b)
    return s1, s2, s3, s4
st = stages(s[a], s[b])
assert st == (sp.Rational(150, 1000), sp.Rational(280, 1000), sp.Rational(300, 1000), sp.Rational(380, 1000)), st
# H+ 消費量からも
assert sp.Rational(3, 10) - (7*sp.Rational(2, 100) + 2*sp.Rational(5, 1000)) == st[0]
# 割合を変えても終点一定（スズ総量 0.0400）
for frac in [0, sp.Rational(1, 4), sp.Rational(1, 2), 1]:
    aa = sp.Rational(4, 100)*frac; bb = sp.Rational(4, 100) - aa
    assert stages(aa, bb)[3] == sp.Rational(38, 100)
pH = -math.log10(2.5e-5 / 99)
assert abs(pH - 6.60) < 0.01 and round(pH, 1) == 6.6
print('INO-AMP-014 個別検証OK', st, round(pH, 3))
