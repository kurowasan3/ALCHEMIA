# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import sqrt
ITEM = '/home/claude/chem2/out/items/THE-ION-043/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

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
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and ' + ' in b:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

def alpha(C, Ka):
    return (-Ka + sqrt(Ka*Ka + 4*Ka*C)) / (2*C)
Ka = 2.0e-4 * 0.5**2 / (1 - 0.5); print('(1)', Ka); assert abs(Ka - 1.0e-4) < 1e-12
a2 = sqrt(Ka/0.10); print('(2) 近似', a2, '厳密', alpha(0.10, Ka)); assert round(a2, 3) == 0.032
a3 = alpha(1.0e-4, Ka); print('(3)', a3, (-1 + 2.24)/2); assert round(a3, 2) == 0.62 and round((-1+2.24)/2, 2) == 0.62
assert abs(sqrt(Ka/1.0e-4) - 1.0) < 1e-12
# (4) 平行移動
KaB = 1.0e-6
print('(4) HB α(2e-6)=', alpha(2.0e-6, KaB))
assert abs(alpha(2.0e-6, KaB) - 0.5) < 1e-9
for C in (1e-1, 1e-3, 1e-5):
    assert abs(alpha(C, Ka) - alpha(C/100, KaB)) < 1e-9
# 図の SVG の規格
svg = open('/home/claude/chem2/out/items/THE-ION-043/problem.svg', encoding='utf-8').read()
assert 'viewBox' in svg and '<script' not in svg and 'foreignObject' not in svg and 'href' not in svg
print('全検証OK:', ITEM)
