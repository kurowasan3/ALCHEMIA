# 検証スクリプト THE-BND-040  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-BND-040/item.json'
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
# 黒鉛の AB 積層を座標で作り，表面層で真下に原子のない原子の割合・最短距離・面密度を数える
d = 0.142; s3 = math.sqrt(3)
a1 = (s3*d, 0.0); a2 = (s3*d/2, 1.5*d)
def layer(shift, N=30):
    return [(i*a1[0] + j*a2[0] + b[0] + shift[0], j*a2[1] + b[1] + shift[1])
            for i in range(-N, N) for j in range(-N, N) for b in [(0, 0), (0, d)]]
A = layer((0, 0)); B = layer((0, d))          # B が表面（第1層），A が第2層
# 中央付近だけで数える
box = lambda p: abs(p[0]) < 2.0 and abs(p[1]) < 2.0
Bc = [p for p in B if box(p)]
Aset = {(round(x, 6), round(y, 6)) for x, y in A}
def over(p): return (round(p[0], 6), round(p[1], 6)) in Aset
# 表面層の原子の半分が A の原子の真上，残りは A の六角形の中心の真上
bright = [p for p in Bc if not over(p)]
frac = len(bright) / len(Bc); print('明るい原子の割合', frac)
assert abs(frac - 0.5) < 0.02
# 真下に原子がない原子は六角形の中心の真上（A の原子 6 個から距離 d）
for p in bright[:20]:
    assert sum(1 for q in A if abs(math.dist(p, q) - d) < 1e-6) == 6
# 明るい点の最短距離と配位数（三角格子なら 6）
p0 = min(bright, key=lambda p: p[0]**2 + p[1]**2)
dists = sorted(math.dist(p0, q) for q in bright if q != p0)
print('最短距離', dists[0], '個数', sum(1 for x in dists if abs(x - dists[0]) < 1e-6))
assert abs(dists[0] - s3*d) < 1e-9 and sum(1 for x in dists if abs(x - dists[0]) < 1e-6) == 6
assert sig(1.73*0.142, 3) == 0.246
# 面密度
area_b = 0.865 * (1.73*0.142)**2
print('明点1個の面積', area_b, '1nm2あたり', 1/area_b)
assert sig(1/area_b, 2) == 19
area_hex = 1.5*1.73*0.142**2; atoms = 2/area_hex
print('原子 /nm2', atoms, '比', (1/area_b)/atoms)
assert abs((1/area_b)/atoms - 0.5) < 0.01
# 数え上げでも確認（4 nm 四方）
assert abs(len([p for p in B if box(p)]) / 16.0 - atoms) / atoms < 0.05
# (4) AA 積層：すべての原子の真下に原子
Aa = layer((0, 0))
assert all(over(p) for p in [q for q in Aa if box(q)])
check_equations()
print('全検証OK:', ITEM)
