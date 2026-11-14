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
Vs = [0, 2, 4, 6, 8, 10, 12, 14, 16]
Ts = [20.000, 20.110, 20.212, 20.306, 20.394, 20.476, 20.461, 20.446, 20.433]
# 表の値はモデル（ΔH1=-12.0 kJ/mol，中和点 10.0 mL）から生成したものと一致するか
for V, T in zip(Vs, Ts):
    n = min(V, 10)*1e-3
    assert abs(20 + 12000*n/((50 + V)*4.2) - T) < 0.0006, (V, T)
# (1) 最高温度の位置
Veq = Vs[Ts.index(max(Ts))]; c = 1.00*Veq/50.0; print('(1)', Veq, c); assert Veq == 10 and rnd(c, 3) == 0.2
# 前半の各点で ΔH がほぼ一定（反応が滴下量に比例）
for V, T in zip(Vs[1:6], Ts[1:6]):
    q = (50 + V)*4.2*(T - 20); dh = q/(V*1e-3)
    assert abs(dh - 12000) < 120, (V, dh)
# (2)
q = 60.0*4.2*(Ts[5] - 20); dH1 = -q/1e-2/1000; print('(2)', dH1); assert rnd(dH1, 3) == -12.0
# (3)
dH2 = dH1 + 56.5; print('(3)', dH2); assert rnd(dH2, 3) == 44.5
dT_HCl = 1e-2*56500/(60*4.2); print('  HCl dT', dT_HCl); assert rnd(dT_HCl, 2) == 2.2
# (4)
Kb = 1e-14/6.0e-10; C = 1e-2/0.06
x = symbols('x', positive=True)
xs = [s for s in solve(Eq(x**2/(C - x), Kb), x) if s > 0][0]
frac = xs/C; print('(4) 加水分解の割合（厳密）', float(frac)); assert rnd(frac, 2) == 0.010
true = 120/(1e-2*(1 - float(frac)))/1000; print('  真の |ΔH1|', true); assert rnd(true, 3) == 12.1
print('ALL OK')
