# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-PER-010/item.json'
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

d = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/elements.json', encoding='utf-8'))
els = [e for v in d.values() if isinstance(v, list) for e in v if isinstance(e, dict) and 'sym' in e]
E = {e['sym']: e for e in els}
pos = {'ア': (3, 13), 'イ': (4, 17), 'ウ': (2, 15), 'エ': (5, 17), 'オ': (4, 11)}
sym = {}
for k, (p, g) in pos.items():
    s = [e['sym'] for e in els if e.get('period') == p and e.get('group') == g]
    assert len(s) == 1; sym[k] = s[0]
print(sym)
assert sym == {'ア': 'Al', 'イ': 'Br', 'ウ': 'N', 'エ': 'I', 'オ': 'Cu'}
states = {k: E[s]['state'] for k, s in sym.items()}; print(states)
assert states['ウ'] == '気体' and states['イ'] == '液体'
assert all(states[k] == '固体' for k in 'アエオ')
ie = {k: E[s]['ie1'] for k, s in sym.items()}; print(ie)
assert max(ie, key=ie.get) == 'ウ'
liq = [e['sym'] for e in els if e.get('state') == '液体' and e.get('period', 9) <= 5]
print('液体(第5周期まで):', liq); assert liq == ['Br']
print('Cu cat:', E['Cu'].get('cat'))

print('全検証OK:', ITEM)
