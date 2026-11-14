# 検証スクリプト THE-BND-039  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-BND-039/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
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
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
import numpy as np, math
r = 1.0; d = 2*r
def stack(seq, nrep=3, n=6):
    """最密層を seq の順に積んだ球の中心座標（層内は 1 辺 2r の三角格子）"""
    a1 = np.array([d,0,0]); a2 = np.array([d/2, d*math.sqrt(3)/2, 0])
    shift = {'A':np.zeros(3), 'B':(a1+a2)/3, 'C':2*(a1+a2)/3}
    h = d*math.sqrt(2/3); pts = []
    for k, L in enumerate(seq*nrep):
        for i in range(-n,n+1):
            for j in range(-n,n+1):
                pts.append(i*a1+j*a2+shift[L]+np.array([0,0,k*h]))
    return np.array(pts), h
for seq in ('ABC','AB'):
    P, h = stack(seq)
    k0 = len(seq)   # 中央付近の層の原点の球を選ぶ
    c = [p for p in P if abs(p[2]-k0*h) < 1e-9 and np.linalg.norm(p[:2]) < 1e-9][0]
    dist = np.linalg.norm(P-c, axis=1)
    nn = P[(abs(dist-2*r) < 1e-6)]
    dz = np.round((nn[:,2]-c[2])/h).astype(int)
    cnt = {z: int((dz==z).sum()) for z in (-1,0,1)}
    print(seq, 'NN', len(nn), cnt)
    assert len(nn) == 12 and cnt == {-1:3, 0:6, 1:3}
    above = [p for p in P if np.linalg.norm(p[:2]-c[:2]) < 1e-9 and p[2] > c[2]+1e-9]
    up = min(p[2]-c[2] for p in above); print(' 真上の最も近い球', up)
    if seq == 'ABC': assert abs(up - 2*math.sqrt(6)*r) < 1e-9
    else: assert abs(up - 4*math.sqrt(6)/3*r) < 1e-9
assert abs(h - 2*math.sqrt(6)/3*r) < 1e-12
# 面心立方格子（立方格子で生成）で体対角線＝2√6 r，2 番目に近い距離と個数
a = 2*math.sqrt(2)*r
basis = np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
F = np.array([(np.array([i,j,k])+b)*a for i in range(-3,4) for j in range(-3,4) for k in range(-3,4) for b in basis])
dist = np.linalg.norm(F, axis=1); ds = sorted(set(np.round(dist,6)))
print('fcc 距離', ds[:4])
assert abs(ds[1]-2*r) < 1e-6 and abs(ds[2]-2*math.sqrt(2)*r) < 1e-6
assert int((abs(dist-ds[2]) < 1e-6).sum()) == 6
assert abs(math.sqrt(3)*a - 2*math.sqrt(6)*r) < 1e-12
# (111) 面の間隔 = a/√3 = h
assert abs(a/math.sqrt(3) - h) < 1e-12
assert abs(2*2.45 - 4.9) < 1e-9

print('ALL OK')
