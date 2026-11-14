# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/out/items/ORG-HYD-039/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# 小問の数と答えの数の一致
nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 0), (nprob, nans)
print('小問数', nprob, '答え', nans)

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def products(n, k):  # 三重結合が k 番目と k+1 番目の炭素の間
    out = set()
    sites = [k + 1] if k == 0 else [k, k + 1]      # 末端は内側の炭素（Markovnikov）
    if n == 2: sites = [0]                           # アセチレン→アセトアルデヒド
    for c in sites:
        atoms = ['C'] * n
        smi = ''.join('C(=O)' if i == c else 'C' for i in range(n))
        out.add(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))
    return out
iodo = Chem.MolFromSmarts('[CH3]C(=O)[#6,#1]')
methylcarbonyl = Chem.MolFromSmarts('[CH3][CX3]=O')
res = {}
for n in range(2, 9):
    for k in range(0, n // 2 if n % 2 == 0 else n // 2 + 1):
        alk = Chem.MolFromSmiles('C' * k + 'C#C' + 'C' * (n - k - 2))
        assert CalcMolFormula(alk) == f'C{n}H{2*n-2}'
        P = products(n, k)
        pos = [Chem.MolFromSmiles(p).HasSubstructMatch(methylcarbonyl) for p in P]
        res[(n, k)] = (P, pos)
        print(n, k, P, pos)
# C6
assert res[(6,0)][0] == {'CCCCC(C)=O'} and len(res[(6,1)][0]) == 2 and res[(6,2)][0] == {'CCCC(=O)CC'}
# C7 で1種類はk=0（1-ヘプチン）だけ，陽性
one7 = [k for (n,k),(P,pos) in res.items() if n == 7 and len(P) == 1]; assert one7 == [0] and all(res[(7,0)][1])
# 1種類かつ陰性（ケトンのみ）
ans = [(n,k) for (n,k),(P,pos) in res.items() if len(P) == 1 and not any(pos) and all('=O)' in p or p.count('C(=O)') for p in P)]
print(ans); assert ans == [(6,2),(8,3)]

print('ALL OK')
