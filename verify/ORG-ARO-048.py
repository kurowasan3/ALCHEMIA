# 検証スクリプト ORG-ARO-048  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ARO-048/item.json'
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
import sympy as sp
from chemlib import parse_formula
AW = {'H': 1.0, 'C': 12, 'N': 14, 'O': 16}
def M(f):
    d = parse_formula(f)[0]; return sum(AW[k]*v for k, v in d.items())
assert (M('C20H15NO3'), M('C7H6O2'), M('C13H11NO2'), M('C6H7NO')) == (317, 122, 213, 109)
nA0 = 317/M('C20H15NO3'); nBz = 134.2/122; nD = 149.1/213
a, b = sp.symbols('a b')
s = sp.solve([a + nD + b - nA0, nD + 2*b - nBz], [a, b])
assert abs(s[a] - 0.100) < 1e-9 and abs(s[b] - 0.200) < 1e-9
assert round(float(s[a])*317, 1) == 31.7 and round(float(s[b])*109, 1) == 21.8
assert round((nD + float(s[b]))/nA0*100) == 90 and round(float(s[b])/nA0*100) == 20
print('収支計算 OK')

print('ALL OK')
