# 検証スクリプト ORG-CBX-024  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/ORG-CBX-024/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, formula_mass
from sympy import symbols, solve, Rational as R
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
                eqs.add(b)
    return sorted(eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
n_items = len(re.findall(r'\\item\b', it['problem'])); n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans == 4, (n_items, n_ans)
eqs = all_equations(); print(eqs)
assert len(eqs) == 2
for e in eqs:
    assert check_equation(e) is None, (e, check_equation(e))
M = {'H': 1.0, 'C': 12, 'O': 16, 'I': 127}
Mst = formula_mass('C57H110O6', M); assert Mst == 890
assert formula_mass('C57H106O6', M) == 886  # ステアリン酸2・リノール酸1
nCC = R(1, 2) * R(1, 10) * R(20, 1000); assert nCC == R(1, 1000)
iv = float(nCC * 254 * 100 / R(443, 1000)); print('ヨウ素価', iv); assert sig(iv, 3) == 57.3
k = symbols('k'); sol = solve(R(443, 1000) / (890 - 4*k) * 2*k - nCC, k); print('k', sol); assert sol == [1]
nIBr = R(1, 2) * R(1, 10) * R(30, 1000); assert nIBr == R(3, 2000)
left = nIBr - nCC; assert left == R(1, 2000); assert left * 2 * 10 * 1000 == 10  # 試料の滴定量 10.0 mL
nCC_big = R(443, 100) / 886 * 2; print('4.43 g の C=C', float(nCC_big)); assert nCC_big == R(1, 100)
assert nCC_big > nIBr
iv_bad = float(nIBr * 254 * 100 / R(443, 100)); print('頭打ちヨウ素価', iv_bad); assert sig(iv_bad, 3) == 8.60
print('全検証OK:', ITEM)
