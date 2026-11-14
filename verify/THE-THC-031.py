# 検証スクリプト THE-THC-031  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-THC-031/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
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
def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', b): continue   # 一般式（R・n）は数えない
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n += 1
    print('反応式', n, '本 OK')
    return n
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    assert it['field'] in FIELDS, it['field']
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)

# ---- 個別の検証 ----
import itertools, math
import sympy as sp
check_structure()
check_all_equations()
# (1) オイラーの定理から記号計算
n, p, h = sp.symbols('n p h')
s = sp.solve([sp.Eq(n - 3*n/2 + p + h, 2), sp.Eq(5*p + 6*h, 3*n)], [p, h])
assert s[p] == 12 and sp.simplify(s[h] - (n - 20)/2) == 0
assert s[h].subs(n, 60) == 20 and s[h].subs(n, 70) == 25 and s[h].subs(n, 78) == 29
# 切頂二十面体（C60）の座標をつくり，辺・5員環・6員環を実際に数える
phi = (1 + 5 ** 0.5) / 2
base = [(0, 1, 3*phi), (1, 2 + phi, 2*phi), (phi, 2, 2*phi + 1)]
V = set()
for b in base:
    for perm in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):      # 偶置換
        q = [b[i] for i in perm]
        for sg in itertools.product((1, -1), repeat=3):
            V.add(tuple(round(sg[i] * q[i], 6) for i in range(3)))
V = sorted(V); assert len(V) == 60
d = lambda a, b: math.dist(a, b)
E = [(i, j) for i, j in itertools.combinations(range(60), 2) if abs(d(V[i], V[j]) - 2) < 1e-6]
assert len(E) == 90
adj = {i: set() for i in range(60)}
for i, j in E: adj[i].add(j); adj[j].add(i)
assert all(len(adj[i]) == 3 for i in adj)
def cycles(k):
    found = set()
    def dfs(path):
        if len(path) == k:
            if path[0] in adj[path[-1]]:
                c = path[:]
                m = c.index(min(c)); c = c[m:] + c[:m]
                if c[1] > c[-1]: c = [c[0]] + c[1:][::-1]
                found.add(tuple(c))
            return
        for nb in adj[path[-1]]:
            if nb not in path: dfs(path + [nb])
    for v in range(60): dfs([v])
    return found
c5, c6 = cycles(5), cycles(6)
print('C60: 辺', len(E), ' 5員環', len(c5), ' 6員環', len(c6)); assert (len(c5), len(c6)) == (12, 20)
# (2)(3)(4)
dHf = 60*(-394) - (-25980); assert dHf == 2340 and sp.Rational(dHf, 60) == 39
eps = sp.Rational(dHf, 12); assert eps == 195
c70 = 70*(-394) - dHf; assert c70 == -29920
per70 = sp.Rational(2340, 70); assert sig(float(per70), 3) == 33.4
nn = sp.solve(sp.Eq((394*n + 2340)/(12*n), sp.Rational(3533, 100)), n)[0]
print('(4) n =', float(nn)); assert round(float(nn)) == 78 and abs(float(nn) - 78) < 0.2
assert abs((394*78 + 2340)/(12*78) - 35.33) < 0.005
# 文献値：ΔfH C60(cr) ≈ +2.3e3 kJ/mol，C60 の燃焼 ≈ −2.60e4 kJ/mol（モデルの値と同程度）
print('THE-THC-031 OK')
