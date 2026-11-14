# 検証スクリプト THE-BND-041  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-BND-041/item.json'
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
import numpy as np, math
def vec(m, ang):   # 置換基側が負を正の値とする（向き：環の中心→置換基）
    return m*np.array([math.cos(math.radians(ang)), math.sin(math.radians(ang))])
def mu(a, b, th): return float(np.linalg.norm(vec(a,0)+vec(b,th)))
TH = {'o':60,'m':120,'p':180}
d = {k: mu(1.70,1.70,t) for k,t in TH.items()}; print('diCl', d)
assert sig(d['o'],2) == 2.9 and sig(d['m'],2) == 1.7 and abs(d['p']) < 1e-9
assert abs(math.sqrt(3)*1.70 - d['o']) < 1e-9
n = {k: mu(1.70,4.00,t) for k,t in TH.items()}; print('ClNO2', n)
assert [round(n[k],2) for k in 'omp'] == [5.07,3.48,2.30]
assert min(n, key=lambda k: abs(n[k]-2.5)) == 'p'
# 他候補との差が十分大きい（一意性）
assert sorted(abs(n[k]-2.5) for k in n)[1] > 0.9
# (3) メチル基の向き
neg = mu(1.70, 0.40, 180); pos = mu(1.70, -0.40, 180); print(neg, pos)
assert abs(neg-1.30) < 1e-9 and abs(pos-2.10) < 1e-9
assert abs(pos-2.1) < abs(neg-2.1)
t = {k: mu(1.70,-0.40,TH[k]) for k in 'om'}; print('ClCH3', t)
assert sig(t['m'],2) == 1.9 and sig(t['o'],2) == 1.5
# (5)
c = 2.50**2/(2*1.70**2) - 1; th = math.degrees(math.acos(c)); print(c, th)
assert round(th) == 85
# 与えた余弦の値
assert abs(math.cos(math.radians(80))-0.174) < 5e-4 and abs(math.cos(math.radians(85))-0.087) < 5e-4
# RDKit で異性体の構造が確かに異なる（o/m/p の正規化 SMILES が 3 種）
from rdkit import Chem
smi = {'o':'Clc1ccccc1[N+](=O)[O-]','m':'Clc1cccc([N+](=O)[O-])c1','p':'Clc1ccc([N+](=O)[O-])cc1'}
assert len({Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in smi.values()}) == 3

print('ALL OK')
