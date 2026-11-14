# 検証スクリプト THE-ION-068  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/THE-ION-068/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
    assert it['field'] in FIELDS, it['field']
    assert all(f in ('qa','equation','compare','condition','pathway','calc','structure') for f in it['format'].split(','))
    import os
    for k in ('problem_fig','solution_fig'):
        if it[k]: assert os.path.exists(os.path.join(os.path.dirname(ITEM), it[k])), it[k]
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
import math, os
from sympy import symbols, nsolve, log, Float, diff, solve
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import GetFormalCharge
lK1, lK2, lKw = -18.8, 1.5, -14.0
def f(pH):
    A = 10**(-2*pH - lK1); B = 10**(lK2 + lKw + pH)
    return 1/(1 + A + B)
vals = [sig(f(p), 2) for p in (9.4, 13.0, 14.0)]
print('f', vals)
assert vals == [0.5, 0.24, 0.031]
assert f(9.4) >= 0.2 and f(13.0) >= 0.2 and f(14.0) < 0.2
# f = 0.20 の厳密解
p = symbols('p')
expr = 1 + 10**(-2*p + 18.8) + 10**(p - 12.5) - 5
lo = float(nsolve(expr, p, 9)); hi = float(nsolve(expr, p, 13))
print('range', lo, hi)
assert round(lo, 1) == 9.1 and round(hi, 1) == 13.1
assert abs(lo - (18.8 - math.log10(4))/2) < 1e-3 and abs(hi - (12.5 + math.log10(4))) < 1e-3
# 最大
import numpy as np
xs = np.linspace(9, 13, 400001); pm = xs[np.argmax([f(x) for x in xs])]
h3 = (lK1 + lK2 + lKw - math.log10(2))/3
print('max pH', pm, -h3, f(pm))
assert abs(pm + h3) < 1e-3 and round(-h3, 1) == 10.5 and sig(f(pm), 2) == 0.98
assert sig(f(12.0), 2) == 0.76
# 構造と分子式・電荷
H2P = 'O=C1OC(c2ccccc12)(c1ccc(O)cc1)c1ccc(O)cc1'
P2 = '[O-]C(=O)c1ccccc1C(=C1C=CC(=O)C=C1)c1ccc([O-])cc1'
C3 = '[O-]C(=O)c1ccccc1C(O)(c1ccc([O-])cc1)c1ccc([O-])cc1'
fm = {k: (CalcMolFormula(Chem.MolFromSmiles(s)), GetFormalCharge(Chem.MolFromSmiles(s))) for k, s in (('H2P', H2P), ('P2', P2), ('C3', C3))}
print(fm)
assert fm['H2P'] == ('C20H14O4', 0) and fm['P2'] == ('C20H12O4-2', -2) and fm['C3'] == ('C20H13O5-3', -3)
assert check_equation('C20H14O4 -> C20H12O4^2- + 2H+') is None
assert check_equation('C20H12O4^2- + OH- -> C20H13O5^3-') is None
# 中心炭素の混成（P2 は sp2，C3 は sp3）
def center(smi):
    m = Chem.MolFromSmiles(smi)
    for a in m.GetAtoms():
        if a.GetSymbol() == 'C' and not a.IsInRing() and sum(1 for n in a.GetNeighbors() if n.IsInRing()) >= 3:
            return str(a.GetHybridization())
print(center(P2), center(C3))
assert center(P2) == 'SP2' and center(C3) == 'SP3'
s = open(os.path.join(os.path.dirname(ITEM), 'problem.svg'), encoding='utf-8').read()
assert 'viewBox' in s and 'script' not in s

print('ALL OK')
