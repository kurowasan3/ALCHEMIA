# THE-SOL-078 検証：キレート抽出の pH 窓
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
from sympy import Rational as Q, symbols, solve, log, N
check_item_equations('/home/claude/chemx/todai/outh/items/THE-SOL-078/item.json')
pH = symbols('pH')
logK = {'In': (Q('0.9'), 3), 'Cu': (Q('1.8'), 2), 'Zn': (Q('-3.0'), 2)}
def logD(m, lhq=-1): k, n = logK[m]; return k + n*lhq + n*pH
print('(2)', solve(logD('In'), pH)); assert solve(logD('In'), pH)[0] == Q(7, 10)
lo = solve(logD('In') - 2, pH)[0]; hi = solve(logD('Zn') + 2, pH)[0]
print('(3)', N(lo), N(hi)); assert round(float(lo), 2) == 1.37 and hi == Q(3, 2)
lo_exact = solve(logD('In') - log(99, 10), pH)[0]; print('  log99 厳密', N(lo_exact))
cu = solve(logD('Cu') - 2, pH)[0]; inn = solve(logD('In') + 2, pH)[0]
print('(4)', N(cu), N(inn)); assert cu > inn
u = symbols('u')
ucu = solve(Q('1.8') + 2*u - 2, u)[0]; uin = solve(Q('0.9') + 3*u + 2, u)[0]
print('u ranges', N(ucu), N(uin)); assert ucu > uin
print('OK')
