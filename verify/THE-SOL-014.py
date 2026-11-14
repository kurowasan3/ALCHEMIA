# 検証スクリプト（g13 作成）：python3 -I このファイル（作業ディレクトリのモジュールに依存しない）
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, Symbol, solve, Eq, nsimplify
ID = __file__.replace('\\', '/').split('/')[-1][:-3]
ITEM = json.load(open(f'/home/claude/chem/out/items/{ID}/item.json', encoding='utf-8'))
ARR = re.compile(r'<=>|->')


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
            if len(ARR.findall(b)) != 1 or '\\boxed' in b or re.search(r'(?<![A-Za-z])[Rn](?![a-z])', b):
                continue
            r = check_equation(b)
            assert r is None, (key, b, r)
            n += 1
    print(f'反応式 {n} 本：原子数・電荷のつり合い OK')


def sig(x, n):
    """有効数字 n 桁に丸めた値"""
    return float(f'{float(x):.{n}g}')


def has(key, *frags):
    for f in frags:
        assert f in ITEM[key], (key, f)


def n_answers():
    a = re.findall(r'(?:^|(?<=[\s　]))\((\d+)\)', ITEM['answer'])
    k = ITEM['problem'].count('\\item')
    assert a == [str(i) for i in range(1, k + 1)], (a, k)
    print(f'小問 {k} 個と答えの番号が一致')

# ---------------------------------------------------------------- 本体
check_all_equations()
n_answers()
a, h = Q(142), Q(180); H = a + h; assert H == 322
fa, fw = a / H, h / H
x = Symbol('x')
# (1)(i) 加熱で無水物
assert 50 - 44 == 6
# (1)(ii) 十水和物
r = solve(Eq((50 - fa * x) / (100 - fw * x), Q(28, 100)), x)[0]; print('1ii', float(r)); assert sig(r, 2) == 77
# (2)(i)(ii)
assert 44 - 44 * Q(50, 100) == 22
y = solve(Eq((44 - fa * x) / (100 - fw * x), Q(28, 100)), x)[0]; print('2ii', float(y), 'Na2SO4', float(y * fa))
assert sig(y, 2) == 56 and y * fa > 22
# (3) 十水和物 32.2 g = 0.100 mol
n = Q('32.2') / H; s_ = n * a; w_ = n * h; print('3', float(s_), float(w_), float(s_ - w_ * Q(50, 100)))
assert s_ == Q('14.2') and w_ == 18 and sig(s_ - w_ / 2, 2) == 5.2
# 実測：Na2SO4 溶解度 25 ℃ 28.1，32.4 ℃（転移点）49.7，80 ℃ 43.7 g/100 g 水．十水和物は 32.38 ℃ で不調和融解
assert check_equation('Na2SO4.10H2O -> Na2SO4 + 10H2O') is None
has('answer', '$6.0\\,\\mathrm g$', '$77\\,\\mathrm g$', '$22\\,\\mathrm g$', '$56\\,\\mathrm g$', '$5.2\\,\\mathrm g$')
print('THE-SOL-014 OK')
