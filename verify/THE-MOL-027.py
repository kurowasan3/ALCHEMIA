# 検証スクリプト THE-MOL-027  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-MOL-027/item.json'
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
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
import sympy as sp
E = sp.symbols('E', positive=True)
sol = sp.solve(sp.Eq(4.75*(E+48), 6.00*(E+35.5)), E); Ev = float(sol[0]); print('X/n', Ev)
assert sig(Ev,3) == 12.0
eq = 4.75/(Ev+35.5); assert abs(eq-0.100) < 1e-12
eqc = 4.20/(Ev+30.0); V = eqc/2*22.4; ox = eqc*(Ev+8.0); print(V, ox)
assert sig(V,3) == 1.12 and sig(ox,3) == 2.00
assert abs(eq*Ev-1.20) < 1e-9 and abs(eq*143.5-14.35) < 1e-9
# 価数 1〜3 で実在の金属か（解説の注）
cands = {n: Ev*n for n in (1,2,3)}; print(cands)   # 12(該当なし), 24(Mg), 36(該当なし)
assert check_equation('MgCO3 -> MgO + CO2') is None
# Mg の実際の値での確認：MgCl2 4.75 g → MgSO4
from rdkit.Chem import Descriptors
print('MgSO4/MgCl2 (24,35.5,32,16):', (24+96)/(24+71)*4.75)
assert abs((24+96)/(24+71)*4.75 - 6.00) < 1e-9
assert abs(4.20/(24+60)*(24+16) - 2.00) < 1e-9

print('ALL OK')
