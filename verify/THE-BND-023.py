# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-023/item.json'
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

rho = {'ice': 0.917, 'water0': 0.99984, 'water4': 0.99997, 'water10': 0.99970, 'CO2_s': 1.56, 'CO2_l': 1.18}
assert rho['ice'] < rho['water0']                     # 氷は浮く
assert rho['water4'] > rho['water0'] and rho['water4'] > rho['water10']   # 4 ℃で最大
# 加圧：体積の小さい（密度の大きい）相へ
assert ('融解' if rho['water0'] > rho['ice'] else '凝固') == '融解'
assert ('融解' if rho['CO2_l'] > rho['CO2_s'] else '凝固') == '凝固'
print('OK')

print('全検証OK:', ITEM)
