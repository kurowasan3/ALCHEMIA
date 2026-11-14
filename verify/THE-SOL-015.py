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
N10, O10, N40, O40 = Q('8.2e-4'), Q('1.7e-3'), Q('5.2e-4'), Q('1.0e-3')
n = N10 * Q(4, 5); o = O10 * Q(1, 5); fr = o / (n + o) * 100
print('N2', float(n), 'O2', float(o), 'O2%', float(fr)); assert sig(n, 2) == 6.6e-4 and sig(o, 2) == 3.4e-4 and sig(fr, 2) == 34
rel = (n - N40 * Q(4, 5)) + (o - O40 * Q(1, 5)); print('released', float(rel)); assert sig(rel, 2) == 3.8e-4
P = Symbol('P')
PN = solve(Eq(N40 * Q(4, 5) * P / Q('1e5'), n), P)[0]; PO = solve(Eq(O40 * Q(1, 5) * P / Q('1e5'), o), P)[0]
print('P_N2', float(PN), 'P_O2', float(PO)); assert PO > PN and sig(PO, 2) == 1.7e5 and sig(PN, 3) == 1.58e5
# 実測（Bunsen 係数から換算，1.013e5 Pa→1.0e5 Pa）：N2 10 ℃ 0.0186，40 ℃ 0.0118；O2 10 ℃ 0.0380，40 ℃ 0.0231 (L/L)
for b, v in [(0.0186, N10), (0.0118, N40), (0.0380, O10), (0.0231, O40)]:
    calc = b / 22.414 / 1.01325; assert abs(calc - float(v)) / calc < 0.05, (b, calc, v)
has('answer', '$6.6\\times10^{-4}\\,\\mathrm{mol}$', '$3.4\\times10^{-4}\\,\\mathrm{mol}$', '$34\\,\\%$', '$3.8\\times10^{-4}\\,\\mathrm{mol}$', '$1.7\\times10^{5}\\,\\mathrm{Pa}$', '酸素')
print('THE-SOL-015 OK')
