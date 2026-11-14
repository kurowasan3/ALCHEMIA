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
ts = [0, 15, 30, 60, 90, 120]; rs = [150.0, 126.6, 110.0, 90.0, 80.0, 75.0]
xa = Rational(70 - 30, 120); print('(1)', xa, (1 - xa)/xa); assert xa == Rational(1, 3) and (1 - xa)/xa == 2
# ずれの半減期を最小二乗で（ln(r-70) の傾き）
import math
xs = ts; ys = [math.log(r - 70) for r in rs]
n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
slope = sum((x - mx)*(y - my) for x, y in zip(xs, ys))/sum((x - mx)**2 for x in xs)
print('slope', slope, 'k1+k2(真の ln)', -slope); assert abs(-slope - math.log(2)/30) < 1e-4
ksum = Rational(693, 1000)/30; k1 = ksum*Rational(2, 3); k2 = ksum/3
print('(2)', float(ksum), float(k1), float(k2)); assert rnd(ksum, 3) == 2.31e-2 and rnd(k1, 3) == 1.54e-2 and rnd(k2, 3) == 7.70e-3
# 速度式を数値積分してデータを再現
a, b, dt = 1.0, 0.0, 0.001; k1f, k2f = math.log(2)/30*2/3, math.log(2)/30/3
t = 0.0; rec = {}
for step in range(int(120/dt) + 1):
    tt = round(step*dt, 3)
    if tt in ts:
        rec[tt] = 150*a + 30*b
    da = (k1f*a - k2f*b)*dt; a -= da; b += da
for t_, r_ in zip(ts, rs):
    assert abs(rec[t_] - r_) < 0.06, (t_, rec[t_], r_)
print('数値積分で表を再現')
# (3)
tt = symbols('tt', positive=True)
th = solve(Eq(Rational(1, 3) + Rational(2, 3)*Rational(1, 2)**(tt/30), Rational(1, 2)), tt)[0]; print('(3)', th); assert th == 60
naive = Rational(693, 1000)/k1; print('  ln2/k1', float(naive)); assert rnd(naive, 2) == 45
# (4)
rb = lambda t_: 70 - 40*Rational(1, 2)**(Rational(t_, 30))
assert rb(30) == 50 and rb(60) == 60
print('ALL OK')
