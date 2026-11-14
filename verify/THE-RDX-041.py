# 検証スクリプト THE-RDX-041  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/THE-RDX-041/item.json'
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
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', b): continue  # 一般式（validate と同じ除外）
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式', n, '本 OK')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure(); check_all_equations()
import sympy as sp
x, y, z = sp.symbols('x y z')
s = sp.solve([x-2.00e-2, y-1.20e-3, 2*x+y+2*z-5.00e-2], [x,y,z])
X, Y, Z = float(s[x]), float(s[y]), float(s[z])
assert sig(Z,3) == 4.40e-3
assert sig(2*X+Y+Z,3) == 4.56e-2 and sig(4*X+5*Y+Z,3) == 9.04e-2
assert sig(X/(X+2*Y)*100,3) == 89.3 and sig(2*X/5.00e-2*100,3) == 80.0
assert sig(2*1.50e-2 + (2.00e-2-1.50e-2),3) == 3.50e-2
sol2 = [(a,b,c,d) for a in range(1,5) for b in range(1,4) for c in range(1,8) for d in range(1,6)
        if check_equation(f'{a}Cu^2+ + {b}HCHO + {c}OH- -> Cu2O + {b}HCOO- + {d}H2O') is None]
print(sol2); assert sol2[0] == (2,1,5,3)
sol4 = [(b,c,d) for b in range(1,4) for c in range(1,8) for d in range(1,6)
        if check_equation(f'Cu^2+ + {b}HCHO + {c}OH- -> Cu + {b}HCOO- + {d}H2O') is None]
print(sol4); assert sol4 == [(1,3,2)]
print('ALL OK')
