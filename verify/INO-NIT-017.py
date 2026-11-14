# 検証スクリプト INO-NIT-017  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/INO-NIT-017/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, formula_mass
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
from fractions import Fraction as F
M = {'H': F(1), 'C': F(12), 'N': F(14), 'O': F(16), 'S': F(32), 'Cl': F(71, 2)}
fs = {'(NH4)2SO4': 2, 'NH4Cl': 1, 'NH4NO3': 2, 'CO(NH2)2': 2}
pct = {}
for f, nN in fs.items():
    m = formula_mass(f, M); pct[f] = 14 * nN / m * 100
    print(f, m, round(float(pct[f]), 2))
order = sorted(pct, key=lambda k: -pct[k])
assert order == ['CO(NH2)2', 'NH4NO3', 'NH4Cl', '(NH4)2SO4'], order
assert [sig(float(pct[k]), 3) for k in ['(NH4)2SO4', 'NH4Cl', 'NH4NO3', 'CO(NH2)2']] == [21.2, 26.2, 35.0, 46.7]
# 残る N（アンモニウムの N だけ失われる）
remain = {'(NH4)2SO4': 0, 'NH4Cl': 0, 'NH4NO3': 100 * F(14) / formula_mass('NH4NO3', M), 'CO(NH2)2': 100 * F(28) / formula_mass('CO(NH2)2', M)}
print({k: float(v) for k, v in remain.items()})
assert sig(float(remain['NH4NO3']), 3) == 17.5 and sig(float(remain['CO(NH2)2']), 3) == 46.7
assert '17.5' in it['answer'] and '46.7' in it['answer']
rx = open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json', encoding='utf-8').read()
assert '2NH4Cl + Ca(OH)2 -> CaCl2 + 2H2O + 2NH3' in rx or '2NH4Cl + Ca(OH)2 -> CaCl2 + 2NH3 + 2H2O' in rx
assert len(all_equations()) == 3

print('ALL OK')
