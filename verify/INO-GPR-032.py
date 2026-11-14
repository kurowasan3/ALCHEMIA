ITEM = '/home/claude/chemx/tokodai/out/items/INO-GPR-032/item.json'
# 検証スクリプト（東工大 g04）  python3 -I この.py
import sys, json, re, math
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
# ---- 個別の検証 ----

from fractions import Fraction as F
rx = {
 'ア':('Cu + 4HNO3 -> Cu(NO3)2 + 2NO2 + 2H2O','NO2'),
 'イ':('3Cu + 8HNO3 -> 3Cu(NO3)2 + 2NO + 4H2O','NO'),
 'ウ':('Cu + 2H2SO4 -> CuSO4 + SO2 + 2H2O','SO2'),
 'エ':('FeS + H2SO4 -> FeSO4 + H2S','H2S'),
 'オ':('2Al + 2NaOH + 6H2O -> 2Na[Al(OH)4] + 3H2','H2'),
 'カ':('CaC2 + 2H2O -> Ca(OH)2 + C2H2','C2H2'),
 'キ':('2NH4Cl + Ca(OH)2 -> CaCl2 + 2NH3 + 2H2O','NH3'),
 'ク':('Mg3N2 + 6H2O -> 3Mg(OH)2 + 2NH3','NH3'),
}
def coef(side, sp):
    for t in side.split(' + '):
        t = t.strip(); m = re.match(r'(\d*)(.*)', t)
        if m.group(2) == sp: return int(m.group(1) or 1)
    return 0
w = {}
for k,(e,g) in rx.items():
    assert check_equation(e) is None, e
    L,R = e.split(' -> ')
    w[k] = F(coef(R,'H2O') - coef(L,'H2O'), coef(R,g))
print(w)
assert w == {'ア':1,'イ':2,'ウ':2,'エ':0,'オ':-2,'カ':-2,'キ':1,'ク':-3}
assert [k for k in w if w[k] < 0] == ['オ','カ','ク']
assert min(w, key=w.get) == 'ク'

print('ALL OK')
