# 検証スクリプト THE-BND-044  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-BND-044/item.json'
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
import itertools
d = 0.142; s3 = math.sqrt(3)
# (1) 六角形あたりの原子数：周期的な網目で数える（単位胞＝原子2個・六角形1個）
a1 = (s3*d, 0.0); a2 = (s3*d/2, 1.5*d)
# (2) 六角形の中心は a1,a2 の三角格子。距離を列挙
cent = [(i*a1[0]+j*a2[0], j*a2[1]) for i in range(-4, 5) for j in range(-4, 5)]
ds = sorted({round(math.hypot(*c), 6) for c in cent if c != (0.0, 0.0)})
print('中心間距離', ds[:3]); assert sig(ds[0], 3) == 0.246 and sig(ds[1], 3) == 0.426
assert ds[0] < 0.40 <= ds[1]
# 最大独立集合（周期境界 L×L のトーラス上の三角格子，最近接を禁止）を厳密に求める
def max_indep(L):
    sites = [(i, j) for i in range(L) for j in range(L)]
    nb = lambda i, j: {((i+di) % L, (j+dj) % L) for di, dj in [(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]}
    best = 0
    # 分枝限定
    order = sites
    def rec(k, chosen, banned):
        nonlocal best
        if len(chosen) + (len(order) - k) <= best: return
        if k == len(order):
            best = max(best, len(chosen)); return
        s = order[k]
        if s not in banned:
            rec(k+1, chosen | {s}, banned | nb(*s) | {s})
        rec(k+1, chosen, banned)
    rec(0, frozenset(), frozenset())
    return best
for L in (3, 6):
    m = max_indep(L); print('L=', L, '最大', m, '/', L*L); assert m * 3 == L*L
# √3×√3 配置の Li–Li 最短距離
occ = [(i*a1[0]+j*a2[0], j*a2[1]) for i in range(-6, 7) for j in range(-6, 7) if (i - j) % 3 == 0]
md = min(math.hypot(p[0]-q[0], p[1]-q[1]) for p, q in itertools.combinations(occ, 2))
print('配置の最短距離', md); assert abs(md - 3*d) < 1e-9
# (3) LiC6
x = 2 * 3; assert x == 6
# (4) 密度
NA = 6.02e23
area_c = 0.5 * 1.5 * 1.73 * 0.142**2
rho_g = 12 / (NA * area_c * 0.335 * 1e-21)
rho_l = (6.9 + 72) / (NA * 6 * area_c * 0.370 * 1e-21)
print('面積/C', area_c, '黒鉛', rho_g, 'LiC6', rho_l)
assert sig(rho_g, 3) == 2.27 and sig(rho_l, 3) == 2.26
inc = (0.370/0.335 - 1) * 100; print('体積増加%', inc); assert sig(inc, 2) == 10
print('質量増加%', (78.9/72 - 1)*100)
# (5) 電気量
Q = 1.00/72 * 9.65e4; print('Q', Q); assert sig(Q, 3) == 1340
assert abs((1/2)/(1/6) - 3) < 1e-12
check_equations()
print('全検証OK:', ITEM)
