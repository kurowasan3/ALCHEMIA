# 検証スクリプト（g03 作成）：python3 -I このファイル
from sympy import *
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ID = __file__.replace('\\', '/').split('/')[-1][:-3]
ITEM = json.load(open(f'/home/claude/chem/out2/items/{ID}/item.json', encoding='utf-8'))
ARR = re.compile(r'<=>|->')
GENERIC = re.compile(r'(?<![A-Za-z])[RnX](?![a-z])')


def ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0:
            return out
        k, d = j + 4, 1
        while d:
            d += {'{': 1, '}': -1}.get(s[k], 0)
            k += 1
        out.append(s[j + 4:k - 1])
        i = k


def check_all_equations():
    """問題文・答え・解説の \\ce{…} のうち矢印1つの完全な式をすべて検査"""
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(ITEM[key]):
            if len(ARR.findall(b)) != 1 or '\\boxed' in b or GENERIC.search(b):
                continue
            r = check_equation(b)
            assert r is None, (key, b, r)
            n += 1
    print(f'反応式 {n} 本：原子数・電荷のつり合い OK')


def check_subitems():
    """小問の数と答えの番号の数が一致するか"""
    k = ITEM['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', ITEM['answer'])
    assert [int(x) for x in nums] == list(range(1, k + 1)), (k, nums)
    print(f'小問 {k} 個と答えの番号が一致')


def near(a, b, rel=0.005):
    a, b = float(a), float(b)
    assert abs(a - b) <= rel * abs(b) + 1e-12, (a, b)


def rnd(x, n):
    """有効数字 n 桁に丸めた値"""
    x = float(x)
    return float(f'{x:.{n}g}')

# ---------------------------------------------------------------- 本体
check_all_equations(); check_subitems()
FM = {'NaCl': 58.5, 'KCl': 74.5, 'NaNO3': 85.0, 'KNO3': 101.0}
ION = {'NaCl': ('Na', 'Cl'), 'KCl': ('K', 'Cl'), 'NaNO3': ('Na', 'NO3'), 'KNO3': ('K', 'NO3')}
SOL = {100: {'NaCl': 39.0, 'KCl': 56.0, 'NaNO3': 180, 'KNO3': 246},
       10: {'NaCl': 35.7, 'KCl': 31.0, 'NaNO3': 80.0, 'KNO3': 22.0}}


def equilibrate(ions, W, T, solid):
    """モデル：少ない方のイオンの量が溶ける量を超えた塩を，超えた分だけ析出（収束まで反復）"""
    ions = dict(ions)
    for _ in range(200):
        moved = False
        for s, (c, an) in ION.items():
            cap = SOL[T][s]/FM[s]*W/100
            ex = min(ions[c], ions[an]) - cap
            if ex > 1e-12:
                ions[c] -= ex; ions[an] -= ex; solid[s] = solid.get(s, 0) + ex; moved = True
        if not moved:
            return ions
    raise RuntimeError


# (1) 300 g から 0.01 g ずつ蒸発させて析出の順を追う
ions = {'Na': 1.0, 'K': 1.0, 'Cl': 1.0, 'NO3': 1.0}; solid = {}; W = 300.0; first = {}
while W > 30:
    ions = equilibrate(ions, W, 100, solid)
    for s, v in solid.items():
        if v > 1e-9 and s not in first:
            first[s] = W
    W = round(W - 0.01, 2)
print('(1) 析出開始の W', first)
assert list(first)[0] == 'NaCl' and abs(first['NaCl'] - 150) < 0.02
assert 'KCl' not in first and 'NaNO3' not in first
assert rnd(first['KNO3'], 3) == 41.1 or abs(first['KNO3'] - 101/246*100) < 0.02
print('KNO3 境界（厳密）', 101/246*100); assert rnd(101/246*100, 3) == 41.1
# (2)
ions = {'Na': 1.0, 'K': 1.0, 'Cl': 1.0, 'NO3': 1.0}; solid = {}
ions = equilibrate(ions, 100, 100, solid)
m2 = solid['NaCl']*58.5; print('(2)', m2, solid); assert set(solid) == {'NaCl'} and rnd(m2, 3) == 19.5
# (3) ろ液のイオン量のままで 10 ℃ の判定（KNO3 析出前の安全側判定）
risk = []
for sname, (c, an) in ION.items():
    if sname == 'KNO3':
        continue
    cap = SOL[10][sname]/FM[sname]
    if min(ions[c], ions[an]) > cap + 1e-12:
        risk.append(sname)
print('(3) おそれ', risk); assert sorted(risk) == ['KCl', 'NaCl']
# 冷却の順序で結果が変わること（モデル上の不定性）の確認：KCl を先に判定すると KCl が出る
s3 = {}; equilibrate(ions, 100, 10, s3); print('  順序依存の例', s3)
# (4) 両方を防ぐ最小の水
Wn = ions['Cl']/(SOL[10]['NaCl']/FM['NaCl'])*100
Wk = min(ions['K'], ions['Cl'])/(SOL[10]['KCl']/FM['KCl'])*100
W2 = max(Wn, Wk); add = W2 - 100
print('(4) W(NaCl)', Wn, 'W(KCl)', Wk, 'add', add); assert rnd(add, 2) == 60
# W2 では，どの順序で判定しても KNO3 以外は出ない
import itertools
for order in itertools.permutations(ION):
    ION2 = {k: ION[k] for k in order}
    st = {}; ii = dict(ions)
    for _ in range(50):
        for sname, (c, an) in ION2.items():
            cap = SOL[10][sname]/FM[sname]*W2/100
            ex = min(ii[c], ii[an]) - cap
            if ex > 1e-9:
                ii[c] -= ex; ii[an] -= ex; st[sname] = st.get(sname, 0) + ex
    assert set(st) == {'KNO3'}, (order, st)
y = 101 - 22.0*W2/100; print('KNO3', y); assert rnd(y, 3) == 65.8
# 中間の温度での単調性（溶解度は温度が高いほど大きい）を表データで確認
assert all(SOL[100][k] > SOL[10][k] for k in FM)
print('ALL OK')
