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

# --- 数値の検証 ---
from sympy import symbols, solve, Rational, nsimplify, sqrt
pM, pD = symbols('pM pD', positive=True)
PA, Pf, R, T, V = 1.50e4, 3.50e4, 8.31e3, 360, 1.00
pu = Pf/2
s = solve([pM+pD-PA, pM+2*pD-pu], [pM, pD])
print(s)
assert abs(s[pM]-1.25e4) < 1e-6 and abs(s[pD]-2.5e3) < 1e-6
Kp = s[pD]/s[pM]**2; print('Kp', Kp); assert abs(Kp-1.60e-5) < 1e-9
frac = 2*s[pD]/pu; print('frac', frac); assert round(float(frac)*100, 1) == 28.6
n = pu*V/(R*T); print('n', n, 'g', n*46); assert round(n*46, 3) == 0.269
# (5) イ：酢酸は 1 単位→1 分子
assert pu == 1.75e4
# (6) 容積半分
x = symbols('x', positive=True)
sol = [r for r in solve(2*Kp*x**2 + x - 2*pu, x) if r > 0][0]
pD2 = Kp*sol**2; frac2 = 2*pD2/(2*pu)
print('half volume: pM', sol, 'pD', pD2, 'frac', frac2)
assert frac2 > frac
check_item('/home/claude/chemx/tokodai/out/items/THE-EQL-124/item.json')
print('OK')
