# 検証スクリプト（グループ g12）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_all_eqs(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)

def check_parts(it):
    n = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n + 1)], (n, nums)
    print('小問数OK:', n)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

ITEM = '/home/claude/chemk/out/items/THE-EQL-149/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_all_eqs(it); check_parts(it)
from sympy import symbols, binomial, simplify, Rational as R
a, k, G, H0 = symbols('a k G H0', positive=True)
H = [binomial(3, i) * a**i * H0 for i in range(4)]
K = [simplify((H[i] * G / H[i+1]).subs(a, G / k)) for i in range(3)]
print('K =', K); assert [simplify(K[0] - k/3), simplify(K[1] - k), simplify(K[2] - 3*k)] == [0, 0, 0]
theta = simplify(sum(i * H[i] for i in range(4)) / (3 * sum(H)))
print('θ =', theta); assert simplify(theta - a/(1+a)) == 0
fr = [simplify(h / sum(H)).subs(a, 1) for h in H]; print('a=1 の割合', fr)
assert fr == [R(1, 8), R(3, 8), R(3, 8), R(1, 8)]
def judge(K1, K2, K3):
    r1, r2 = K2 / K1, K3 / K2
    if abs(r1 - 3) < 1e-9 and abs(r2 - 3) < 1e-9: return '独立'
    if r1 < 3 and r2 < 3: return '正'
    if r1 > 3 and r2 > 3: return '負'
    return '混合'
X = (R(2, 10**4), R(6, 10**4), R(18, 10**4)); Y = (R(2, 10**4),)*3; Z = (R(1, 10**4), R(1, 10**3), R(1, 10**2))
res = [judge(*X), judge(*Y), judge(*Z)]; print(res); assert res == ['独立', '正', '負']
kX = 3 * X[0]; assert kX == R(6, 10**4) and X[1] == kX and X[2] == 3 * kX
aa = R(12, 10**4) / kX; th = aa / (1 + aa); print('θ =', float(th)); assert sig(float(th), 2) == 0.67
print('ALL OK')
