import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def load(item_id):
    return json.load(open(f'/home/claude/chem/out2/items/{item_id}/item.json', encoding='utf-8'))

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

def check_all_equations(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return len(eqs)

def check_subitems(it):
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    if n_items:
        assert [int(x) for x in nums] == list(range(1, n_items + 1)), (n_items, nums)
    print('小問数OK:', n_items)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def in_text(val_str, it, key='answer'):
    assert val_str in it[key], (val_str, key)
# ---- POL-SUG-019 グルコースの変旋光（独立計算）
from sympy import Rational as R, nsimplify
it = load('POL-SUG-019'); check_all_equations(it); check_subitems(it)
a, b, eq, gl = 112.0, 18.7, 52.7, 158.0
y = (a - eq) / (a - b)                 # beta fraction
print('beta', y, 'alpha', 1 - y, 'K', y / (1 - y))
assert sig(100*y, 3) == 63.6 and sig(100*(1-y), 3) == 36.4 and sig(y/(1-y), 3) == 1.74
for s_ in ['$36.4\\,\\%$', '$63.6\\,\\%$', '$K=1.74$']: in_text(s_, it)
# (3)
w = (gl - 135.0) / (gl - a); mg = 2.00 * w
fin = w * eq + (1 - w) * gl
M_glc = 6*12 + 12*1.0 + 6*16; M_cu2o = 2*63.5 + 16
cu = mg / M_glc * M_cu2o
print('(3)', mg, fin, cu)
assert sig(mg, 3) == 1.00 and round(fin) == 105 and sig(cu, 3) == 0.794
# (4)
w2 = (gl - 135.0) / (gl - b); mb = 2.00 * w2
fin2 = w2 * eq + (1 - w2) * gl; cu2 = mb / M_glc * M_cu2o
print('(4)', w2, mb, fin2, cu2)
assert 0 < w2 < 1 and sig(mb, 3) == 0.330 and round(fin2) == 141 and sig(cu2, 3) == 0.262
assert fin < 135 < fin2   # 減少 / 増加
# 図鑑データ: グルコースの還元性・メチルグリコシド
d = open('/home/claude/kurowasan3/alchemia/content/chem/organic.json', encoding='utf-8').read()
assert 'グルコース' in d
print('POL-SUG-019 OK')
