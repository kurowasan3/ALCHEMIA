# 検証スクリプト ORG-ARO-080  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/outh/items/ORG-ARO-080/item.json'
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
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


from rdkit import Chem
cs = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
def ring_positions(m, cN):
    dm = Chem.GetDistanceMatrix(m)
    ring = [r for r in m.GetRingInfo().AtomRings() if cN in r][0]
    return ring, dm
def cN_of(m): return [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='C' and any(n.GetSymbol()=='N' for n in a.GetNeighbors())][0]
def add_br(m, idxs):
    rw = Chem.RWMol(m)
    for i in idxs:
        k = rw.AddAtom(Chem.Atom(35)); rw.AddBond(i, k, Chem.BondType.SINGLE)
    mm = rw.GetMol(); Chem.SanitizeMol(mm); return mm
def r1(smi):
    m = Chem.MolFromSmiles(smi); c = cN_of(m); ring, dm = ring_positions(m, c)
    return add_br(m, [i for i in ring if dm[c][i] in (1,3) and m.GetAtomWithIdx(i).GetTotalNumHs()==1])
def r2(smi):
    m = Chem.MolFromSmiles(smi); c = cN_of(m); ring, dm = ring_positions(m, c)
    free = lambda d: [i for i in ring if dm[c][i]==d and m.GetAtomWithIdx(i).GetTotalNumHs()==1]
    t = free(3) or free(1)
    outs = {Chem.MolToSmiles(add_br(m, [i])) for i in t}
    assert len(outs) == 1, outs
    return add_br(m, [t[0]])
def r3(m):
    rw = Chem.RWMol(m); n = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='N'][0]
    rw.RemoveAtom(n); mm = rw.GetMol(); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
ANI, PT, OT, MT = 'Nc1ccccc1', 'Nc1ccc(C)cc1', 'Nc1ccccc1C', 'Nc1cccc(C)c1'
a1 = r3(r1(ANI)); assert a1 == cs('Brc1cc(Br)cc(Br)c1')
a2 = r3(r1(PT)); a3 = r3(r1(OT)); print(a2, a3)
assert a2 == a3 == cs('Cc1cc(Br)cc(Br)c1')
a4 = r3(r2(OT)); print(a4); assert a4 == cs('Cc1cccc(Br)c1')
assert r3(r2(PT)) == cs('Cc1cccc(Br)c1')   # 原題の経路と同じ生成物になることの確認
tri = [s for s in (OT, MT, PT) if r3(r1(s)).count('Br') == 3]; print(tri)
assert tri == [MT] and r3(r1(MT)) == cs('Cc1c(Br)cc(Br)cc1Br')
print('OK')
