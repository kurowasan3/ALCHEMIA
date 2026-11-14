# 検証スクリプト INO-HAL-015  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/INO-HAL-015/item.json'
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
M = {'H': F(1), 'O': F(16), 'Cl': F(71, 2), 'Ca': F(40)}
m1 = formula_mass('Ca(ClO)2.2H2O', M); m2 = formula_mass('CaCl(ClO).H2O', M)
assert m1 == 179 and m2 == 145, (m1, m2)
n1 = F(358, 100) / m1 * 2; n2 = F(358, 100) / m2
v1 = n1 * F(224, 10); v2 = n2 * F(224, 10)
print('V', float(v1), float(v2), 'ratio', float(n1 / n2))
assert sig(float(v1), 2) == 0.90 and sig(float(v2), 2) == 0.55 and sig(float(n1 / n2), 2) == 1.6
for e in ['Cl2 + 2NaOH -> NaCl + NaClO + H2O', 'NaClO + 2HCl -> NaCl + Cl2 + H2O', 'ClO- + Cl- + 2H+ -> Cl2 + H2O',
          'Ca(ClO)2.2H2O + 4HCl -> CaCl2 + 2Cl2 + 4H2O', 'CaCl(ClO).H2O + 2HCl -> CaCl2 + Cl2 + 2H2O']:
    assert e in all_equations() or e.replace(' ', '') in [x.replace(' ', '') for x in all_equations()], e
# 図鑑との照合
rx = open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json', encoding='utf-8').read()
assert 'CaCl(ClO).H2O + 2HCl -> CaCl2 + 2H2O + Cl2' in rx and 'Cl2 + 2NaOH -> NaCl + NaClO + H2O' in rx
assert '0.90' in it['answer'] and '0.55' in it['answer'] and '1.6' in it['answer']

print('ALL OK')
