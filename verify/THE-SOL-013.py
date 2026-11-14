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
d = Q('1.38'); w = Q('0.630'); M = Q('63.0')
c = 1000 * d * w / M; m = (1000 * w / M) / ((1000 - 1000 * w) / 1000)
print('mol/L', float(c), 'mol/kg', float(m)); assert sig(c, 3) == 13.8 and sig(m, 3) == 27.0
md = Q('0.100') / ((1000 - Q('0.100') * M) / 1000); print('dilute mol/kg', float(md)); assert sig(md, 3) == 0.101
mass = Q('1.12') * 100; acid = mass * Q('0.200'); A = acid / w; VA = A / d; Vw = mass - A
print('A g', float(A), 'A mL', float(VA), 'water mL', float(Vw), 'sum', float(VA + Vw))
assert sig(VA, 3) == 25.8 and sig(Vw, 3) == 76.4 and VA + Vw > 100
# 実測密度（20 ℃）：HNO3 63 % ≈ 1.38，20 % ≈ 1.115 g/mL → 混合で体積は約 2 % 縮む（現実的）
has('answer', '$13.8\\,\\mathrm{mol/L}$', '$27.0\\,\\mathrm{mol/kg}$', '$0.101\\,\\mathrm{mol/kg}$', '$25.8\\,\\mathrm{mL}$', '$76.4\\,\\mathrm{mL}$')
print('THE-SOL-013 OK')
