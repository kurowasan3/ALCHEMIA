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

from sympy import symbols, log, Rational, series, exp, N
import math
c0, V, nO2, dt = 0.800, 0.0200, 4.00e-3, 5.00
c1 = c0 - 2*nO2/V; assert abs(c1-0.400) < 1e-12
kapp = ((c0-c1)/dt)/((c0+c1)/2); k = math.log(2)/dt
print('kapp', kapp, 'k', k, 'ratio', kapp/k, '<c>', (c0-c1)/(k*dt))
assert round(kapp, 3) == 0.133 and round(k, 3) == 0.139 and round(kapp/k, 3) == 0.962
# 問題で与えた値でも
assert round((2/3)/0.693, 3) == 0.962
ratio = lambda r: (2*(1-r)/(1+r))/(-math.log(r))
vals = {r: ratio(r) for r in (0.9, 0.7, 0.5, 0.25)}; print(vals)
ok = [r for r, v in vals.items() if 1 - v <= 0.02]; assert ok == [0.9, 0.7]
# つねに 1 未満
assert all(ratio(i/1000) < 1 for i in range(1, 1000))
x = symbols('x', positive=True)
print('series:', series(2*(1-exp(-x))/(1+exp(-x))/x, x, 0, 4))
# (5)
step = lambda c, d: c*(1-kapp*d/2)/(1+kapp*d/2)
a = step(step(c0, 5), 5); b = step(c0, 10); true = c0*math.exp(-k*10)
print('method A', a, 'method B', b, 'true', true)
assert round(a, 3) == 0.200 and round(b, 3) == 0.160 and abs(true-0.2) < 1e-12
check_item('/home/claude/chemx/tokodai/outh/items/THE-EQL-130/item.json')
print('OK')
