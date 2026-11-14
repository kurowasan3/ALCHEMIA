# 検証スクリプト THE-MOL-028  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-MOL-028/item.json'
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
from math import gcd
ratio = {'I':17.0,'II':33.0,'III':68.75,'IV':62.0}
wX = {'I':0.912,'II':0.939,'III':0.225,'IV':1.0}
M = {k: 2.0*v for k,v in ratio.items()}
mX = {k: M[k]*wX[k] for k in M}
print(M, mX)
assert [sig(M[k],3) for k in 'I II III IV'.split()] == [34.0,66.0,138,124]
assert [sig(mX[k],3) for k in 'I II III IV'.split()] == [31.0,62.0,30.9,124]
# 公約数 31
r = [round(mX[k]/31) for k in mX]; assert r == [1,2,1,4]
assert all(abs(mX[k]/31 - round(mX[k]/31)) < 0.01 for k in mX)
# デュロン・プティ
for m,val in [(1,23.9),(2,11.9),(3,8.0)]:
    assert abs(0.77*31/m - val) < 0.06, (m, 0.77*31/m)
# 分子式（H=1.0, Cl=35.5, X=31）
assert abs(M['I'] - (31+3*1.0)) < 1e-9 and abs(M['II'] - (62+4)) < 1e-9 and abs(M['III'] - (31+3*35.5)) < 1e-9 and M['IV'] == 124
# 実在物質（P）との整合：RDKit で質量パーセントを確認
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
for smi, n, w in [('P', 1, 0.912), ('PP', 2, 0.939), ('ClP(Cl)Cl', 1, 0.225)]:
    mw = MolWt(Chem.MolFromSmiles(smi)); frac = n*30.974/mw
    print(smi, round(mw,2), round(frac,4)); assert abs(frac - w) < 0.002
# (5) 単原子水素の仮定
M1 = {k: 1.0*v for k,v in ratio.items()}
H_I = M1['I']*(1-0.912); print('H in I', H_I); assert abs(H_I-1.5) < 0.01
mX1 = {k: M1[k]*wX[k] for k in M1}; print(mX1)
assert abs(0.77*15.5 - 11.9) < 0.05
# 【A】の式
assert check_equation('H2 + Cl2 -> 2HCl') is None

print('ALL OK')
