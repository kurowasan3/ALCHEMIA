CODE = 'THE-GAS-029'
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

from fractions import Fraction as Fr
T  = [80, 85, 90, 95, 100, 105, 110, 120, 130]
PW = [Fr(x).limit_denominator(100) for x in ['4.7','5.8','7.0','8.5','10.0','12.1','14.3','19.9','27.0']]
PX = [Fr(x).limit_denominator(100) for x in ['1.9','2.4','3.0','3.6','4.3','5.1','6.0','8.0','10.0']]
# 表が問題文と一致
for t, w, x in zip(T, PW, PX):
    assert f' {t} ' in J['problem'] or f'{t} &' in J['problem']
S = [w + x for w, x in zip(PW, PX)]
assert S[2] == 10 and S[1] < 10 and all(s > 10 for s in S[3:])
assert all(PW[i] < PW[i+1] and PX[i] < PX[i+1] for i in range(8))
assert PX[-1] == 10 and PW[4] == 10   # 純液体の沸点
mass_ratio = Fr(3) * Fr(225, 2) / (Fr(7) * 18)
assert round(float(mass_ratio), 1) == 2.7
cb_out = Fr(1) * Fr(3, 7)
assert cb_out < Fr(1, 2) and Fr(1, 2) * Fr(7, 3) > 1
left = (Fr(1, 2) - cb_out) * Fr(225, 2)
assert round(float(left), 1) == 8.0, float(left)
R = [Fr(2, 3) * w + Fr(1, 3) * x for w, x in zip(PW, PX)]
assert R[5] < 10 < R[6], (float(R[5]), float(R[6]))
assert R[4] < 10
y_lo, y_hi = PX[5] / 30, PX[6] / 30
assert 0.16 < y_lo < 0.18 and 0.19 < y_hi < 0.21 and y_hi < Fr(3, 10)
print('THE-GAS-029 個別検証OK', float(mass_ratio), float(left), float(R[5]), float(R[6]), float(y_lo), float(y_hi))
