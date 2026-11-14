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

ITEM = '/home/claude/chemk/out/items/THE-ACB-043/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_all_eqs(it); check_parts(it)
pKa = {'P': 0.99, 'Q': -4.5, 'R': -6.6, 'S': -10.1}
def ratio(pk, H0): return 10 ** (pk - H0)
def color(r): return '酸性色' if r >= 10 else ('塩基性色' if r <= 0.1 else '中間色')
r1 = ratio(pKa['P'], 2.00); print('(1)', r1, sig(r1, 2), color(r1))
assert sig(r1, 2) == 0.098 and color(r1) == '塩基性色'
c2 = {k: color(ratio(v, -5.0)) for k, v in pKa.items()}; print('(2)', c2)
assert c2 == {'P': '酸性色', 'Q': '中間色', 'R': '塩基性色', 'S': '塩基性色'}
def bracket(changed):
    up = min(pKa[k] for k in changed)            # H0 <= 変色した中で最小の pKa
    low = max(v for k, v in pKa.items() if k not in changed)  # H0 > 変色しない中で最大の pKa
    return low, up
X = bracket({'P', 'Q'}); Y = bracket({'P', 'Q', 'R'}); print('(3)', X, Y)
assert X == (-6.6, -4.5) and Y == (-10.1, -6.6)
r4 = ratio(pKa['R'], -8.0); fr = r4 / (1 + r4) * 100; print('(4)', fr)
assert round(fr) == 96
print('ALL OK')
