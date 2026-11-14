# THE-THC-047 検証：イオン半径モデルの溶解熱
# --- 共通：item.json の \ce{} 反応式のつり合い検査（単独で動く） ---
import json, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
def _ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k, d = j + 4, 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k
def check_item_equations(path):
    it = json.load(open(path, encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for ce in _ce_blocks(it[key]):
            arrows = re.findall(r'<=>|->|<-', ce)
            if len(arrows) != 1 or '\\boxed' in ce or '\\text' in ce or '\\mathrm' in ce:
                continue
            r = check_equation(ce)
            assert r is None, (key, ce, r)
            n += 1
    print('反応式', n, '本 OK')
    return it
from sympy import symbols, simplify, Rational as Q, N
check_item_equations('/home/claude/chemx/todai/outh/items/THE-THC-047/item.json')
b, rm, rx = symbols('beta r_M r_X', positive=True)
Qexpr = -4*b/(rm + rx) + b*(1/rm + 1/rx)
assert simplify(Qexpr - b*(rm - rx)**2/(rm*rx*(rm + rx))) == 0
r = {'Li': Q('0.09'), 'K': Q('0.15'), 'Cs': Q('0.17'), 'F': Q('0.12'), 'Cl': Q('0.17'), 'Br': Q('0.18'), 'I': Q('0.21')}
def q(m, x): return float(Qexpr.subs({b: 56, rm: r[m], rx: r[x]}))
vals = {m + x: q(m, x) for m in ['Li', 'Cs', 'K'] for x in ['F', 'Cl', 'Br', 'I']}
for k, v in vals.items(): print(k, round(v, 2))
assert round(vals['LiF']) == 22 and round(vals['LiI'], -1) == 140 and round(vals['CsF']) == 24 and round(vals['CsI'], 1) == 6.6
assert vals['LiI'] > vals['LiF'] and vals['CsF'] > vals['CsI']
assert min(['KF', 'KCl', 'KBr', 'KI'], key=lambda k: vals[k]) == 'KCl'
# 元問題の NaF の値（Q_ion=-923, Q_Na+Q_F=930, r=0.12）から α/β
alpha = 923*Q('0.24'); beta = 930/(2/Q('0.12')); print('alpha/beta', N(alpha/beta))
print('OK')
