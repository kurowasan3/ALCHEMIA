# 検証スクリプト ORG-HYD-031（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it

it = check_item('/home/claude/chem5/out/items/ORG-HYD-031/item.json')

from rdkit import Chem
from itertools import combinations
def chlorinated(smi, k):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); Hs = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'H']; out = set()
    for c in combinations(Hs, k):
        rw = Chem.RWMol(m)
        for i in c: rw.GetAtomWithIdx(i).SetAtomicNum(17)
        mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(Chem.RemoveHs(mm)))
    return out
cnt = [len(chlorinated('CC', k)) for k in range(1, 7)]
print(cnt); assert cnt == [1, 2, 2, 2, 1, 1] and sum(cnt) == 9
pct = {k: round(35.5 * k / (24 + 1.0 * (6 - k) + 35.5 * k) * 100, 1) for k in range(1, 7)}
print(pct); assert [k for k, v in pct.items() if v == 71.7] == [2]
assert check_equation('CH2=CHCl + HCl -> CH3CHCl2') is None

print('ALL OK')
