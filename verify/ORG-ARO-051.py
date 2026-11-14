# 検証スクリプト ORG-ARO-051  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-ARO-051/item.json'
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
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b) and (' + ' in b or re.search(r'(^|\s)\d+\s*[A-Z(\[]', b)):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert len(nums) == n_items or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items, '答え数', len(nums))
check_structure()
check_all_equations()

# ---- 個別の検証 ----

import sympy as sp
from rdkit import Chem
R = sp.Integer(24)
fo = R*6*sp.Rational(58, 100)/2; fm = R*6*sp.Rational(4, 100)/2; fp = R*6*sp.Rational(38, 100)
assert (fo, fm, fp) == (sp.Rational(4176, 100), sp.Rational(288, 100), sp.Rational(5472, 100))
assert sig(float(fo), 2) == 42 and sig(float(fm), 2) == 2.9 and sig(float(fp), 2) == 55
# 位置ごとの関係を環の位置から求める
def rel(p, q):
    d = min((p - q) % 6, (q - p) % 6)
    return {1: 'o', 2: 'm', 3: 'p'}[d]
F = {'o': fo, 'm': fm, 'p': fp}
def predict(me):
    contrib = {}
    for pos in range(1, 7):
        if pos in me: continue
        f = 1
        for q in me: f *= F[rel(pos, q)]
        # 生成物の同一性は RDKit で判定
        mol = Chem.MolFromSmiles(ring(me, pos))
        key = Chem.MolToSmiles(mol)
        contrib[key] = contrib.get(key, 0) + f
    return contrib
def ring(me, pos):
    parts = []
    for i in range(1, 7):
        sub = '(C)' if i in me else ('([N+](=O)[O-])' if i == pos else '')
        parts.append(sub)
    return 'c1' + parts[0] + ''.join('c' + p for p in parts[1:5]) + 'c1' + parts[5]
res = {}
for name, me in {'o': (1, 2), 'm': (1, 3), 'p': (1, 4)}.items():
    c = predict(me)
    tot = sum(c.values())
    res[name] = (tot / 6, sorted((float(v / tot) for v in c.values()), reverse=True))
print({k: (float(v[0]), v[1]) for k, v in res.items()})
assert sig(float(res['p'][0]), 2) == 80 and len(res['p'][1]) == 1
assert sig(float(res['m'][0]), 2) == 1100 and len(res['m'][1]) == 3
mfr = res['m'][1]
assert sig(mfr[0]*100, 2) == 72 and sig(mfr[1]*100, 2) == 28 and sig(mfr[2]*100, 2) == 0.13
assert sig(float(res['o'][0]), 2) == 93 and len(res['o'][1]) == 2
assert sig(res['o'][1][0]*100, 2) == 57
assert abs(res['o'][1][0] / res['o'][1][1] - 4/3) < 0.05

print('ALL OK')
