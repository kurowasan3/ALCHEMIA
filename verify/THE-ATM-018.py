# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem4/out/items/THE-ATM-018/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
def count(atoms, charge):   # atoms: [(元素, 質量数)]
    p = sum(Z[e] for e, A in atoms); n = sum(A - Z[e] for e, A in atoms)
    return p, n, p - charge
sp = {
 'N2': ([('N',14)]*2, 0), 'CO': ([('C',12),('O',16)], 0), 'CN-': ([('C',12),('N',14)], -1),
 'NO+': ([('N',14),('O',16)], 1), 'O2': ([('O',16)]*2, 0), 'OH-': ([('O',16),('H',1)], -1),
 'C2^2-': ([('C',12)]*2, -2), 'D2O': ([('H',2),('H',2),('O',16)], 0), 'H2O18': ([('H',1),('H',1),('O',18)], 0)}
res = {k: count(*v) for k, v in sp.items()}
for k, v in res.items(): print(k, '陽子,中性子,電子 =', v)
same_e = sorted(k for k in ['CO','CN-','NO+','O2','OH-','C2^2-'] if res[k][2] == res['N2'][2])
assert same_e == sorted(['CO','CN-','NO+','C2^2-']), same_e
assert [k for k in same_e if res[k] == res['N2']] == ['CO']
assert res['D2O'] == res['H2O18'] == (10, 10, 10)
print('OK')
