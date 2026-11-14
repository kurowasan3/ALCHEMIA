# 検証スクリプト INO-ALK-018  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/INO-ALK-018/item.json'
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
M = {'H': 1.008, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'Na': 22.99, 'Cl': 35.45, 'K': 39.10}
sol = {'NaCl': (35.9, 58.5), 'NH4Cl': (37.2, 53.5), 'NaHCO3': (9.6, 84.0), 'NH4HCO3': (21.7, 79.0), 'KCl': (34.2, 74.6), 'KHCO3': (33.7, 100.1)}
for f, (s, fw) in sol.items():
    assert abs(formula_mass(f, M) - fw) < 0.15, (f, formula_mass(f, M))
mol = {f: s / fw for f, (s, fw) in sol.items()}
print({k: round(v, 3) for k, v in mol.items()})
na = ['NaCl', 'NH4Cl', 'NaHCO3', 'NH4HCO3']; k = ['KCl', 'NH4Cl', 'KHCO3', 'NH4HCO3']
assert min(na, key=mol.get) == 'NaHCO3' and min(k, key=mol.get) == 'NH4HCO3'
assert [sig(mol[x], 2) for x in na] == [0.61, 0.70, 0.11, 0.27]
assert [sig(mol[x], 2) for x in k] == [0.46, 0.70, 0.34, 0.27]
rx = open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json', encoding='utf-8').read()
assert 'NaCl + H2O + NH3 + CO2 -> NaHCO3 + NH4Cl' in rx
assert len(all_equations()) == 3

print('ALL OK')
