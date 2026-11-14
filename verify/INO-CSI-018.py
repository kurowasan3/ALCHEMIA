# 検証スクリプト INO-CSI-018   python3 -I この.py
ITEM = '/home/claude/chemk/out/items/INO-CSI-018/item.json'
import sys, json, re
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
def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def check_structure():
    n_items = len(re.findall(r'\\item\b', re.sub(r'\\begin\{itemize\}.*?\\end\{itemize\}', '', it['problem'], flags=re.S)))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    nums = [int(x) for x in re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, n_items + 1)), nums
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
check_all_equations()
# ---- 個別の検証（SymPy で物質収支）----
import sympy as sp
Si_, H, Cl = sp.Rational(281, 10), sp.Rational(1), sp.Rational(355, 10)
nSi = sp.Rational(1000) * sp.Rational(98, 100) / Si_
nT = nSi * sp.Rational(9, 10); nC = nSi - nT
mass = nT * Si_ / 1000
print('Si', float(nSi), 'SiHCl3', float(nT), 'SiCl4', float(nC), '析出 kg', float(mass))
assert sig(float(mass), 2) == 0.88
hcl_in = 3 * nT + 4 * nC; hcl_back = 3 * nT; make_up = hcl_in - hcl_back
h2_out = nT + 2 * nC; h2_use = nT; h2_ex = h2_out - h2_use
print('HCl 補給', float(make_up), 'H2 余り', float(h2_ex))
assert sig(float(make_up), 2) == 14 and sig(float(h2_ex), 2) == 7.0
# 正味の変化が Si + 4HCl -> SiCl4 + 2H2 の倍数になっていること
assert make_up == 4 * nC and h2_ex == 2 * nC
# 沸点による分離の判定（SiHCl3 32 ℃）
bp = {'BCl3': 13, 'PCl3': 76}
assert bp['BCl3'] < 32 < bp['PCl3']
print('すべて OK')
