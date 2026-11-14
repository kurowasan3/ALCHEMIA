# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-PER-011/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

T = {'W': [738, 1451, 7733, 10543, 13630, 18020], 'X': [801, 2427, 3660, 25026, 32827],
     'Y': [419, 3052, 4420, 5877, 7975, 9590], 'Z': [1012, 1907, 2914, 4964, 6274, 21268]}
nist = {'Mg': [737.7, 1450.7, 7732.7, 10542.5, 13630, 18020], 'B': [800.6, 2427.1, 3659.7, 25025.8, 32826.7],
        'K': [418.8, 3052, 4420, 5877, 7975, 9590], 'P': [1011.8, 1907, 2914.1, 4963.6, 6273.9, 21267]}
val = {}
for k, v in T.items():
    ratios = [v[i + 1] / v[i] for i in range(len(v) - 1)]
    val[k] = ratios.index(max(ratios)) + 1
print('価電子', val); assert val == {'W': 2, 'X': 3, 'Y': 1, 'Z': 5}
ref = {'Li': 520, 'Na': 496, 'Be': 900, 'Ca': 590, 'Al': 578, 'N': 1402, 'As': 944}
# 同族で周期が大きいほど E1 小
assert ref['Be'] > T['W'][0] > ref['Ca']      # Mg（第3周期）
assert T['X'][0] > ref['Al']                   # Al より上 → B
assert T['Y'][0] < ref['Na'] < ref['Li']       # Na より下 → K（Rb は第5周期で範囲外）
assert ref['N'] > T['Z'][0] > ref['As']        # P
for (k, s) in zip('WXYZ', ['Mg', 'B', 'K', 'P']):
    assert all(abs(a - b) <= 1 for a, b in zip(T[k], nist[s])), k
print('NIST 値と一致')
d = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/elements.json', encoding='utf-8'))
els = {e['sym']: e for v in d.values() if isinstance(v, list) for e in v if isinstance(e, dict) and 'sym' in e}
for s in list(ref) + ['Mg', 'B', 'K', 'P']:
    print(s, els[s]['ie1'], ref.get(s, ''))
    if s in ref: assert abs(els[s]['ie1'] - ref[s]) <= 2

print('全検証OK:', ITEM)
