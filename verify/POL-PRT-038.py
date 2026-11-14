# 検証スクリプト POL-PRT-038  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/POL-PRT-038/item.json'
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
        assert r is None, (e, r)
        print('式OK:', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
from itertools import combinations, combinations_with_replacement as cwr
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
check_structure(); check_all_equations()
AA = {'Gly':(2,5,1,2,0),'Ala':(3,7,1,2,0),'Ser':(3,7,1,3,0),'Cys':(3,7,1,2,1),'Val':(5,11,1,2,0),'Leu':(6,13,1,2,0),
      'Met':(5,11,1,2,1),'Phe':(9,11,1,2,0),'Tyr':(9,11,1,3,0),'Asp':(4,7,1,4,0),'Glu':(5,9,1,4,0),'Lys':(6,14,2,2,0)}
W = (12,1,14,16,32)
def pep(c):
    t = [sum(AA[a][i] for a in c) for i in range(5)]; n = len(c)
    t[1] -= 2*(n-1); t[3] -= n-1; return tuple(t)
target = (17,25,3,5,1)
mw = sum(a*b for a,b in zip(target, W)); assert mw == 383
pct = [round(100*a*b/mw, 1) for a,b in zip(target, W)]; print('質量%', pct)
assert pct == [53.3, 6.5, 11.0, 20.9, 8.4]
# 組成式：S=1 で割った比
r = [x/w for x,w in zip(pct, W)]; print('比', [round(x/r[4],2) for x in r])
assert 2*mw > 500
tri = sorted(tuple(sorted(c)) for c in combinations(AA, 3) if pep(c) == target)
print('異なる3種', tri); assert tri == [('Ala','Met','Tyr'), ('Cys','Tyr','Val'), ('Met','Phe','Ser')]
anyn = sorted({tuple(sorted(c)) for n in range(2,6) for c in cwr(AA, n) if pep(c) == target}); assert anyn == tri
# 和の分子式
s = tuple(sum(AA[a][i] for a in ('Cys','Tyr','Val')) for i in range(5)); assert s == (17,29,3,7,1)
# (3) キサントプロテイン：全候補が Tyr/Phe を含む．ヨウ素：Cys を含むものだけ
assert all({'Tyr','Phe'} & set(c) for c in tri)
assert [c for c in tri if 'Cys' in c] == [('Cys','Tyr','Val')]
assert 2*mw - 2 == 764
q = Chem.MolFromSmiles('N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CS)C(=O)N[C@@H](C(C)C)C(=O)O')
print(CalcMolFormula(q)); assert CalcMolFormula(q) == 'C17H25N3O5S'
# 一般式の反応 2R-SH + I2 -> R-S-S-R + 2HI：R を Q として原子数確認
assert check_equation('2C17H25N3O5S + I2 -> C34H48N6O10S2 + 2HI') is None
print('POL-PRT-038 OK')
