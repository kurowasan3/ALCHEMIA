# 検証スクリプト POL-RES-019  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/POL-RES-019/item.json'
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
from fractions import Fraction as F
m1 = F(940,10) + F(360,10) - 18*F(12,10)
assert m1 == F(1084,10) and sig(float(m1), 3) == 108
assert F(3,2) == F(15,10)
mol = 1 - F(8,10); mass = 94 + 30*F(8,10) - 18*F(8,10)
assert mol == F(2,10) and mass == F(1036,10)
assert sig(float(mass/mol), 2) == 520 and 1/mol == 5
# 環をつくらない結合の数え方をシミュレーションで確認
import random
random.seed(1)
N = 1000; parent = list(range(N)); sites = [3]*N
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
bonds = 0
while bonds < 800:
    a, b = random.randrange(N), random.randrange(N)
    if a == b or sites[a] == 0 or sites[b] == 0 or find(a) == find(b): continue
    parent[find(a)] = find(b); sites[a] -= 1; sites[b] -= 1; bonds += 1
assert len({find(i) for i in range(N)}) == N - 800
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
assert CalcMolFormula(Chem.MolFromSmiles('Oc1ccccc1Cc1ccccc1O')) == 'C13H12O2'  # 2C6H6O + CH2O - H2O
assert check_equation('2C6H6O + CH2O -> C13H12O2 + H2O') is None
print('フェノール樹脂の量的関係 OK')

print('ALL OK')
