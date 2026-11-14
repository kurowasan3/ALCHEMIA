# 検証スクリプト THE-MOL-053   python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-MOL-053/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
import sympy as sp
for e in ['NH3BH3 -> NH2BH2 + H2', 'NH3BH3 -> NHBH + 2H2', '3NH3BH3 -> B3N3H6 + 6H2']:
    assert check_equation(e) is None, e
M = {'AB': 14+10.8+6*1.0, 'a': 14+10.8+4*1.0, 'b': 14+10.8+2*1.0, 'bz': 3*(14+10.8)+6*1.0}
print(M); assert M == {'AB': 30.8, 'a': 28.8, 'b': 26.8, 'bz': 80.4} or all(abs(M[k]-v) < 1e-9 for k, v in {'AB': 30.8, 'a': 28.8, 'b': 26.8, 'bz': 80.4}.items())
n0 = 1.000/M['AB']; H = 1.31/22.4; m = 0.753
a, b, c = sp.symbols('a b c')
s = sp.solve([a+b+c-n0, a+2*b+2*c-H, M['a']*a+M['b']*b-m], [a, b, c])
fr = {k: float(v/n0)*100 for k, v in s.items()}; print('(3)', fr)
assert [sig(fr[a],2), sig(fr[b],2), sig(fr[c],2)] == [20, 65, 15]
# H2 の測定だけでは b,c が決まらない：係数行列の階数
A = sp.Matrix([[1,1,1],[1,2,2]]); assert A.rank() == 2 and A[:, 1] == A[:, 2]
mbz = float(s[c]/3*M['bz']); print('(4)', mbz, m + H*2.0 + mbz); assert sig(mbz,3) == 0.130
assert abs(m + H*2.0 + mbz - 1.000) < 1e-9
r = (1.000-m)/2.0/H; print('(5)', r); assert sig(r,2) == 2.1
# 設定の妥当性：真の割合 20/65/15 % から作ったデータが与えた数値に丸まる
n = 1/30.8; Ht = n*(0.20+2*0.65+2*0.15)*22.4; mt = n*(0.20*28.8+0.65*26.8)
print('元データ', Ht, mt); assert round(Ht,2) == 1.31 and round(mt,3) == 0.753

print('検証OK')
