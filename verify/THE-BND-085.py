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

from sympy import symbols, sqrt, pi, simplify, diff, Rational, factor
rho, r = symbols('rho r', positive=True)
vol = 4*Rational(4,3)*pi*r**3*(1+rho**3)
f1 = simplify(vol/(2*r*(1+rho))**3)
f2 = simplify(vol/(2*sqrt(2)*r)**3)
assert simplify(f1 - Rational(2,3)*pi*(1+rho**3)/(1+rho)**3) == 0
assert simplify(f2 - pi*(1+rho**3)/(3*sqrt(2))) == 0
rc = sqrt(2)-1
# 接触条件: a/√2 >= 2r with a=2r(1+rho) <=> rho >= √2-1
assert simplify((sqrt(2)*(1+rc)) - 2) == 0
assert simplify(f1.subs(rho, rc) - f2.subs(rho, rc)) == 0
g = (1+rho**3)/(1+rho)**3
assert simplify(diff(g, rho) - 3*(rho-1)/(1+rho)**3) == 0
fmax = simplify(f2.subs(rho, rc)); print('fmax', fmax, float(fmax))
assert simplify(fmax - pi*(5-3*sqrt(2))/3) == 0
from sympy import limit
assert simplify(limit(f2, rho, 0) - pi/(3*sqrt(2))) == 0      # 面心立方の充塡率
assert simplify(f1.subs(rho, 1) - pi/6) == 0
# 与えた近似値で
P, s2, s3 = 3.14, 1.414, 1.732
fN = lambda x: (2*P/3)*(1+x**3)/(1+x)**3; fA = lambda x: P*(1+x**3)/(3*s2); fCs = lambda x: (s3*P/2)*(1+x**3)/(1+x)**3
vals = dict(fmax=P*(5-3*s2)/3, LiI=fA(0.076/0.220), KCl=fN(0.138/0.181), lim0=P/(3*s2))
print(vals)
assert round(vals['fmax'], 3) == 0.793 and round(vals['LiI']*100, 1) == 77.1 and round(vals['KCl']*100, 1) == 55.2 and round(vals['lim0'], 3) == 0.740
assert 0.076/0.220 < s2-1 < 0.138/0.181
check_item('/home/claude/chemx/tokodai/outh/items/THE-BND-085/item.json')
print('OK')
