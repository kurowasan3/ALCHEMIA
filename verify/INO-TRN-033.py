# 検証スクリプト INO-TRN-033  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/INO-TRN-033/item.json'
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
NA, a = 6.0e23, 3.62e-8
rho = 4*63.5/(NA*a**3); assert sig(rho,3) == 8.92
t = 0.180/(sig(rho,3)*40); assert sig(t*1e4,2) == 5.0
n = 0.180/63.5; assert sig(n*22.4e3,2) == 63
# 係数の一意性：a,b,c,d を総当たり
sols = [(A,B,C,D) for A in range(1,5) for B in range(1,9) for C in range(0,3) for D in range(0,6)
        if C == 1 and check_equation(f'Cu^2+ + {A}HCHO + {B}OH- -> Cu + {A}HCOO- + {C}H2 + {D}H2O') is None]
print(sols); assert sols == [(2,4,1,2)]
print('ALL OK')
