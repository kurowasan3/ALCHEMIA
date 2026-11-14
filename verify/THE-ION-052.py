# 検証スクリプト THE-ION-052  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ION-052/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    keys = ['src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig']
    assert all(k in it for k in keys), [k for k in keys if k not in it]
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
    print('小問数', n_items, '答え数', len(nums))
check_structure()
check_all_equations()

# ---- 個別の検証 ----

import numpy as np
from scipy.optimize import least_squares
Ksp, KspBr, b2, b4, Kb, Kw = 1.8e-10, 5.0e-13, 1.6e7, 1e9, 1.8e-5, 1e-14
def solve(V=1.0, HNO3=0, Zn=0, NaOH=0, KBr=0):
    Agt, Clt, Nt, Brt = 0.050, 0.050, 0.30, KBr
    NO3 = 0.10 + HNO3 + 2*Zn; Na = NaOH; K = KBr
    def f(z):
        lh, lN, lAg, nCl, nBr, lZn = z
        h, N, Ag, Zn2 = np.exp([lh, lN, lAg, lZn])
        Cl = Ksp/Ag; Br = KspBr/Ag if Brt > 0 else 0
        oh = Kw/h; NH4 = Kb*N/oh
        AgA = b2*Ag*N**2; ZnA = b4*Zn2*N**4 if Zn > 0 else 0
        Zn2 = Zn2 if Zn > 0 else 0
        r = [ (nCl + (nBr if Brt>0 else 0) + V*(Ag+AgA) - Agt)/Agt,
              (nCl + V*Cl - Clt)/Clt,
              (V*(N+NH4+2*AgA+4*ZnA) - Nt)/Nt,
              (h + Na/V + K/V + Ag + AgA + NH4 + 2*Zn2 + 2*ZnA - oh - Cl - Br - NO3/V)/0.1 ]
        r.append((nBr + V*Br - Brt)/1e-3 if Brt > 0 else nBr)
        r.append((V*(Zn2+ZnA) - Zn)/Zn if Zn > 0 else lZn)
        return r
    best = None
    for g in [(-21, np.log(0.17), np.log(3e-8), 0.04, 0.0, -20)]:
        s = least_squares(f, g, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if best is None or s.cost < best.cost: best = s
    assert best.cost < 1e-20, best.cost
    return best.x
base = solve()
nAgCl0 = base[3]
print('基準: AgCl(s) =', nAgCl0, 'pH =', -np.log10(np.exp(base[0])))
assert 0 < nAgCl0 < 0.050
res = {'HNO3': solve(HNO3=0.005)[3], 'Zn': solve(Zn=0.002)[3], 'NaOH': solve(NaOH=0.005)[3],
       'KBr': solve(KBr=0.005)[3], 'dil': solve(V=2.0)[3]}
print(res)
assert res['HNO3'] > nAgCl0 and res['Zn'] > nAgCl0
assert res['NaOH'] < nAgCl0 and res['KBr'] < nAgCl0
assert abs(res['dil'] - nAgCl0)/nAgCl0 < 0.01
print('J129 OK')
