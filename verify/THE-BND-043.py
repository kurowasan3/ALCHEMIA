# 検証スクリプト THE-BND-043  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-BND-043/item.json'
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡'.split()
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
def check_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and all(x.strip() for x in re.split(r'->|<=>|<-', b)):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK', e)
    return sorted(eqs)
def sig(x, n):
    return float(f'{x:.{n}g}')
def check_structure():
    assert it['field'] in FIELDS
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    for f in ('problem_fig', 'solution_fig'):
        if it[f]:
            import os
            p = os.path.join(os.path.dirname(ITEM), it[f]); assert os.path.exists(p), p
            t = open(p).read(); assert 'viewBox' in t and '<script' not in t and 'foreignObject' not in t
    print('小問数', n_items, '答え数', n_ans)
check_structure()
# 4 化合物それぞれ「陽–陰接触（a=2(r+ + r-)，r+/r- >= √2-1）」か「陰–陰接触（a=2√2 r-，r+/r- <= √2-1）」
# 測定値は a±0.0005 nm の幅で扱い，線形計画で実行可能な場合を列挙する
from scipy.optimize import linprog
import itertools
S2 = 1.414; k = S2 - 1
data = {('A','X'): 0.509, ('A','Y'): 0.622, ('B','X'): 0.560, ('B','Y'): 0.640}
var = {'A': 0, 'B': 1, 'X': 2, 'Y': 3}
tol = 0.0005
def feasible(pattern, objective=None):
    Aub, bub = [], []
    for (c, an), mode in zip(data, pattern):
        a = data[(c, an)]
        row = [0]*4
        if mode == 'cat':
            row[var[c]] = 2; row[var[an]] = 2
            r2 = [0]*4; r2[var[c]] = -1; r2[var[an]] = k       # r+ >= k r-
            Aub.append(r2); bub.append(0)
        else:
            row[var[an]] = 2*S2
            r2 = [0]*4; r2[var[c]] = 1; r2[var[an]] = -k        # r+ <= k r-
            Aub.append(r2); bub.append(0)
        Aub.append(row); bub.append(a + tol)
        Aub.append([-x for x in row]); bub.append(-(a - tol))
    res = []
    for i in range(4):
        out = []
        for sgn in (1, -1):
            cvec = [0]*4; cvec[i] = sgn
            r = linprog(cvec, A_ub=Aub, b_ub=bub, bounds=[(1e-4, 1)]*4, method='highs')
            if r.status != 0: return None
            out.append(sgn * r.fun)
        res.append(out)
    return res
keys = list(data)
feas = {}
for pat in itertools.product(['cat', 'an'], repeat=4):
    r = feasible(pat)
    if r: feas[pat] = r; print([f'{c}{a}:{m}' for (c, a), m in zip(keys, pat)], [(round(x[0],4), round(x[1],4)) for x in r])
# すべての実行可能な場合で：AX・AY では陰イオンどうしが（少なくとも）接している＝ a ≒ 2√2 r-
for pat, r in feas.items():
    rX = r[2]; rY = r[3]
    assert abs(data[('A','X')] - 2*S2*rX[0]) < 0.0025 and abs(data[('A','X')] - 2*S2*rX[1]) < 0.0025
    assert abs(data[('A','Y')] - 2*S2*rY[0]) < 0.0025 and abs(data[('A','Y')] - 2*S2*rY[1]) < 0.0025
    # B 化合物は陽–陰接触
    assert pat[2] == 'cat' and pat[3] == 'cat'
assert ('cat','cat','cat','cat') not in feas
# (1) 加成性のチェック
s1 = data[('A','X')] + data[('B','Y')]; s2 = data[('A','Y')] + data[('B','X')]
print('AX+BY', s1, 'AY+BX', s2); assert abs(s1 - s2) > 0.02
# (2) 半径
rX = 0.509/(2*S2); rY = 0.622/(2*S2); rB1 = 0.560/2 - rX; rB2 = 0.640/2 - rY
print('rX', rX, 'rY', rY, 'rB', rB1, rB2)
assert sig(rX, 3) == 0.180 and sig(rY, 3) == 0.220 and sig(rB1, 2) == 0.10 and sig(rB2, 2) == 0.10
assert rB1/rX > k and rB2/rY > k and rB1/rY < 0.732
# (3) rA の上限
rAmax = k * rX; print('rA <=', rAmax); assert round(rAmax, 3) == 0.075 or sig(rAmax, 2) == 0.075
rA_all = [r[0] for r in feas.values()]
print('LP の rA 範囲', min(x[0] for x in rA_all), max(x[1] for x in rA_all))
assert max(x[1] for x in rA_all) < 0.0760
# (4) AZ（rZ = 0.140）の格子定数の範囲
rZ = 0.140
amin = 2*S2*rZ; amax = 2*(round(rAmax, 4) + rZ)
print('AZ a 範囲', amin, amax); assert sig(amin, 3) == 0.396 and sig(amax, 3) == 0.429
# rA = rAmax のとき AZ は陽–陰接触（rA/rZ > k）
assert rAmax/rZ > k
# 境界：rA < k rZ = 0.0580 なら a = 0.396（陰–陰接触）
print('境界 rA', k*rZ)

check_equations()
print('全検証OK:', ITEM)
