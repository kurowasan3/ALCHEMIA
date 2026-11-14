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

ITEM = '/home/claude/chemk/outh/items/THE-ION-119/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_all_eqs(it); check_parts(it)
from math import log10
L2, L3 = 0.30, 0.48
lg = {0.50: -L2, 8.0: 3*L2, 0.125: -3*L2, 6.0: L2+L3}
for k, v in lg.items(): assert abs(log10(k) - v) < 0.01, k
used = [('測定1 B1', 0.50), ('S1 B1', 8.0), ('S1 B2', 0.125), ('S2 B2', 6.0), ('S2 B3', 0.125), ('S3 B3', 8.0)]
for name, r in used: assert 0.1 <= r <= 10, name
print('用いた r はすべて 0.1〜10 の範囲内')
for mode in ('given', 'exact'):
    f = (lambda x: lg[x]) if mode == 'given' else log10
    pKa1 = 1.00 + f(0.50)
    H1 = pKa1 - f(8.0); pKa2 = H1 + f(0.125)
    H2 = pKa2 - f(6.0); pKa3 = H2 + f(0.125)
    H3 = pKa3 - f(8.0)
    vals = [round(x, 2) for x in (pKa1, H1, pKa2, H2, pKa3, H3)]
    print(mode, vals)
    if mode == 'given':
        assert vals == [0.70, -0.20, -1.10, -1.88, -2.78, -3.68], vals
        G = pKa1, H1, pKa2, H2, pKa3, H3
pKa1, H1, pKa2, H2, pKa3, H3 = G
# S3 で B1, B2 は範囲外
for pk in (pKa1, pKa2):
    assert abs(pk - H3) > 1
print('S3 での log r (B1, B2):', round(pKa1 - H3, 2), round(pKa2 - H3, 2))
# Γ
logG = -log10(8.0) - H3 if False else -3*L2 - H3
print('log Γ =', logG, 'Γ =', 10**logG); assert sig(10**logG, 2) == 600
# (5)
lo, hi = pKa2 - 1, pKa2 + 1
print('B2 の測定可能な H0:', lo, hi, '必要な [H+] >=', 10**(-hi), 'r=1 で', 10**(-pKa2))
assert sig(10**(-hi), 2) == 1.3 and round(10**(-pKa2)) == 13
r3 = 10**(pKa3 - H1); print('S1 での B3 の r =', r3); assert sig(r3, 2) == 0.0026
# S1, S2 で測った両方の指示薬が範囲内（自己無撞着）
assert 0.1 <= 10**(pKa1 - H1) <= 10 and 0.1 <= 10**(pKa2 - H1) <= 10
assert 0.1 <= 10**(pKa2 - H2) <= 10 and 0.1 <= 10**(pKa3 - H2) <= 10
assert 0.1 <= 10**(pKa3 - H3) <= 10
print('ALL OK')
