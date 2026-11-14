# 検証スクリプト（グループ g02）  python3 -I この.py
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def check_all_equations(it, extra=()):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|…', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]', b): continue
            if not (' + ' in b or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', b)): continue
            e = clean_ce(b)
            r = check_equation(e)
            assert r is None, (b, r)
            print('式OK:', b); n += 1
    for e in extra:
        r = check_equation(e); assert r is None, (e, r); print('式OK(追加):', e); n += 1
    return n

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def check_parts(it):
    k = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, k + 1)], (k, nums)
    print('小問数OK:', k)
ITEM = '/home/claude/chemk/out/items/ORG-CBX-075/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); assert check_all_equations(it) >= 6
from rdkit import Chem
EN = {'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'P': 2.19, 'Cl': 3.16, 'S': 2.58}
def oxnums(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); Chem.Kekulize(m, clearAromaticFlags=True)
    ox = {a.GetIdx(): -a.GetFormalCharge()*0 for a in m.GetAtoms()}
    for a in m.GetAtoms(): ox[a.GetIdx()] = a.GetFormalCharge()
    for b in m.GetBonds():
        i, j = b.GetBeginAtom(), b.GetEndAtom(); o = int(b.GetBondTypeAsDouble())
        ei, ej = EN[i.GetSymbol()], EN[j.GetSymbol()]
        if ei > ej: ox[i.GetIdx()] -= o; ox[j.GetIdx()] += o
        elif ej > ei: ox[j.GetIdx()] -= o; ox[i.GetIdx()] += o
    return m, ox
def of(smi, sym):
    m, ox = oxnums(smi); return sorted(ox[a.GetIdx()] for a in m.GetAtoms() if a.GetSymbol() == sym)
PPh3, OPPh3 = 'P(c1ccccc1)(c1ccccc1)c1ccccc1', 'O=P(c1ccccc1)(c1ccccc1)c1ccccc1'
DEAD, DEADH2 = 'CCOC(=O)N=NC(=O)OCC', 'CCOC(=O)NNC(=O)OCC'
print('P', of(PPh3, 'P'), '->', of(OPPh3, 'P')); assert of(PPh3, 'P') == [3] and of(OPPh3, 'P') == [5]
print('N', of(DEAD, 'N'), '->', of(DEADH2, 'N')); assert of(DEAD, 'N') == [-1, -1] and of(DEADH2, 'N') == [-2, -2]
# 他の元素の酸化数の総和が不変
lhs = ['CCC(=O)O', 'CO', PPh3, DEAD]; rhs = ['CCC(=O)OC', OPPh3, DEADH2]
def tot(smis, sym): return sum(sum(of(s, sym)) for s in smis)
for sym in ['C', 'H', 'O']:
    print(sym, tot(lhs, sym), tot(rhs, sym)); assert tot(lhs, sym) == tot(rhs, sym)
print('全検証OK:', ITEM)
