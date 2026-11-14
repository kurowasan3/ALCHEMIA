# --- 反応式の機械的チェック（item.json 中の \ce{...} で矢印が1つの完全な式） ---
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def ce_blocks(s):
    out = []
    i = 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0:
            return out
        k = j + 4; depth = 1
        while depth:
            if s[k] == '{': depth += 1
            elif s[k] == '}': depth -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k

def check_item(path):
    d = json.load(open(path, encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(d[key]):
            if '\\text' in b or '\\mathrm' in b or '\\boxed' in b:
                continue
            arrows = re.findall(r'<=>|->|<-', b)
            if len(arrows) != 1:
                continue
            e = re.sub(r'\[[^\]]*\]', '', b)          # 矢印ラベルを除く
            e = e.replace('<=>', '->')
            e = re.sub(r'\((g|l|s|aq)\)', '', e)
            r = check_equation(e)
            assert r is None, (b, r)
            n += 1
    print('反応式チェック OK:', n, '本')

from sympy import symbols, solve, simplify
P, P1, V, w, s = symbols('P P1 V w s', positive=True)
sol = solve(P1*V - P*(V-w) - P*s*w, P)[0]
print('P =', sol); assert simplify(sol - P1*V/(V-w+s*w)) == 0
f = lambda s_: 1.00e5*2.00/(2.00-1.00+s_*1.00)
print('N2', f(0.0160), 'frac', 0.016/1.016, ' CO2', f(1.70))
assert round(f(0.016)/1e5, 2) == 1.97 and round(0.016/1.016*100, 2) == 1.57 and round(f(1.70)/1e4, 2) == 7.41
# P>P1 <=> s<1
expr = simplify(sol/P1 - 1)
print('P/P1-1 =', expr)
assert simplify(expr - (1-s)*w/(V-w+s*w)) == 0
check_item('/home/claude/chemx/tokodai/out/items/THE-SOL-090/item.json')
print('OK')
