# 検証スクリプト POL-SUG-038（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    print('小問の数', np_, na)
    return it

it = check_item('/home/claude/chem6/out/items/POL-SUG-038/item.json')
from sympy import Rational as Q
Mn = Q(100,10000) * Q(840,100) / 1000          # 0.0100 mol/L × 8.40 mL
Fe = 5 * Mn; Cu2O = Fe / 2; sugar10 = Cu2O; sugar = sugar10 * 10
mass = sugar * 180; pct = mass / Q(1,2) * 100
print('MnO4-', float(Mn), 'Fe2+', float(Fe), '糖(全体)', float(sugar), 'g', float(mass), '%', float(pct))
assert float(Q(2,5)) == 0.4 and abs(float(pct) - 75.6) < 1e-9
Cr = Fe / 6; V = Cr / Q(1,100) * 1000
print('Cr2O7 2-', float(Cr), 'mL', float(V)); assert abs(float(V) - 7.00) < 1e-9
print('OK')
