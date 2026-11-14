# 検証スクリプト THE-SOL-039  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-SOL-039/item.json'
import sys, json, re
from math import log10, floor
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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
def check_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '->', 'OK' if r is None else r)
        assert r is None, (e, r)
    return sorted(eqs)
def sig(x, n):
    return float(f'{x:.{n}g}')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def close(a, b, rel=1e-9):
    assert abs(a - b) <= rel * max(abs(a), abs(b)), (a, b)
check_structure()
print('反応式', check_equations())
# 文献値（Shannon 半径 nm，水和エンタルピー kJ/mol，イオンの移動度 1e-8 m^2/(V s)，モル伝導率 S cm^2/mol）
r = {'Li': 0.076, 'Na': 0.102, 'K': 0.138, 'Mg': 0.072}
z = {'Li': 1, 'Na': 1, 'K': 1, 'Mg': 2}
dH = {'Li': -520, 'Na': -405, 'K': -321, 'Mg': -1920}
order_field = sorted(r, key=lambda i: -z[i]/r[i])
order_dH = sorted(r, key=lambda i: dH[i])
print('電荷/半径の順', order_field, '水和エンタルピーの順', order_dH)
assert order_field == order_dH == ['Mg', 'Li', 'Na', 'K']
mob = {'Li': 4.01, 'Na': 5.19, 'K': 7.62}
lam = {'LiCl': 115.0, 'NaCl': 126.5, 'KCl': 149.9}
assert sorted(mob, key=mob.get) == ['Li', 'Na', 'K'] and sorted(lam, key=lam.get) == ['LiCl', 'NaCl', 'KCl']
# 問題文の表の値が文献値と同じこと
for ion, v in (('Li', '0.076'), ('Na', '0.102'), ('K', '0.138'), ('Mg', '0.072')):
    assert v in it['problem']
assert '\\ce{Mg^2+} > \\ce{Li+} > \\ce{Na+} > \\ce{K+}' in it['answer']
print('すべて一致')
