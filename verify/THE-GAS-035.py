# 検証スクリプト THE-GAS-035  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-GAS-035/item.json'
import sys, json, re
from math import log10, floor
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
def check_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '->', 'OK' if r is None else r)
        assert r is None, (e, r)
    return sorted(eqs)
def sig(x, n):
    return float(f'{x:.{n}g}')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def close(a, b, rel=1e-9):
    assert abs(a - b) <= rel * max(abs(a), abs(b)), (a, b)
check_structure()
from math import log, exp
print('反応式', check_equations())
R = 8.31
# (1)
dV = (18.0/1.000 - 18.0/0.917)*1e-6
slope = 6.0e3/(273*dV); print('(1) 傾き', slope); assert sig(slope, 2) == -1.3e7
P = 60*9.8/1.0e-4; dT = P/abs(slope); print('(1) 刃の圧力', P, '融点低下', dT); assert sig(dT, 2) == 0.44 and dT < 5
# (2)
Hsub = 6.0 + 45; assert Hsub == 51
sv = 45e3*610/(R*273**2); ss = 51e3*610/(R*273**2); print('(2) 傾き 蒸気圧', sv, '昇華圧', ss)
assert sig(sv, 2) == 44 and sig(ss, 2) == 50 and ss > sv
x = 1/263 - 1/273
ev, es = -45e3/R*x, -51e3/R*x; print('(2) 指数', ev, es); assert round(ev, 3) == -0.754 and round(es, 3) == -0.855
Pl, Ps = 610*exp(ev), 610*exp(es); print('(2) 水', Pl, '氷', Ps); assert sig(Pl, 2) == 290 and sig(Ps, 2) == 260 and Ps < Pl
assert abs(exp(-0.754) - 0.470) < 5e-4 and abs(exp(-0.855) - 0.425) < 5e-4
# (3) 与えた表の ln 値が P と一致するか
data = {260: (1.17e3, 7.063), 270: (2.48e3, 7.816), 290: (8.28e3, 9.021), 310: (2.06e4, 9.931)}
for T, (p, l) in data.items():
    assert abs(log(p) - l) < 3e-3, (T, log(p), l)
bs = (7.816-7.063)/(1/260-1/270); bl = (9.931-9.021)/(1/290-1/310)
Hs, Hv = bs*R, bl*R; Hf = Hs - Hv
print('(3) 昇華熱', Hs, '蒸発熱', Hv, '融解熱', Hf)
assert sig(Hs/1e3, 2) == 44 and sig(Hv/1e3, 2) == 34 and sig(Hf/1e3, 2) == 9.9
a_s = 7.816 + bs/270; a_l = 9.931 + bl/310
invT = (a_s - a_l)/(bs - bl); Tt = 1/invT; lnPt = a_l - bl*invT
print('(3) 三重点', Tt, exp(lnPt), 'lnP', lnPt); assert sig(Tt, 2) == 280 and sig(exp(lnPt), 2) == 5000
assert abs(log(1.013e5) - 11.526) < 1e-3 and abs(exp(8.52) - 5.0e3) < 50
Tb = bl/(a_l - log(1.013e5)); print('(3) 沸点', Tb); assert sig(Tb, 2) == 350
# 解説中の中間値
print('b_s', bs, 'a_s', a_s, 'b_l', bl, 'a_l', a_l)
assert round(bs) == 5286 and round(bl) == 4090 and round(a_s, 3) == 27.394 and round(a_l, 3) == 23.126
assert round(a_s - a_l, 3) == 4.268 and round(1/Tt, 6) == 0.003570 and round(lnPt, 2) == 8.52
print('すべて一致')
