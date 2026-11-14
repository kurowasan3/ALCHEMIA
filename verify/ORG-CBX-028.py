# ORG-CBX-028 検証：グリセリンのエステル（トリ・ジグリセリド）の構造異性体・立体異性体の数
from itertools import product
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
ACY = {'H': 'O', 'a': 'OC(=O)C', 'b': 'OC(=O)CC', 'c': 'OC(=O)CCC'}
def mol(g):
    return Chem.MolFromSmiles(f"C({ACY[g[0]]})C({ACY[g[1]]})C({ACY[g[2]]})")
def counts(assigns):
    flat, stereo, chiral = set(), set(), set()
    for g in assigns:
        m = mol(g)
        flat.add(Chem.MolToSmiles(m, isomericSmiles=False))
        for s in EnumerateStereoisomers(m):
            cs = Chem.MolToSmiles(s)
            stereo.add(cs)
            mir = Chem.MolToSmiles(Chem.MolFromSmiles(cs.replace('@@', '!').replace('@', '@@').replace('!', '@')))
            if mir != cs: chiral.add(cs)
    return len(flat), len(stereo), len(chiral)
acids = ['a', 'b', 'c']
for n in (1, 2, 3):
    A = acids[:n]
    f, s, c = counts(product(A, repeat=3))
    print('tri n=', n, f, s, c)
    assert f == n*n*(n+1)//2 and s == n**3 and c == n**3 - n**2
    # ジグリセリド（OH 1 個）
    dg = [g for g in product(['H'] + A, repeat=3) if g.count('H') == 1]
    f2, s2, c2 = counts(dg)
    print('di  n=', n, f2, s2, c2)
    assert f2 == n*n + n*(n+1)//2 and s2 == 3*n*n
# (5) P, O, S を 1 個ずつ：2 位が S の候補（鏡像異性体を区別）
cand = sorted(g for g in product('POS', repeat=3) if sorted(g) == ['O', 'P', 'S'] and g[1] == 'S')
assert cand == [('O', 'S', 'P'), ('P', 'S', 'O')]
# 位置 1 だけを切る酵素で O が遊離 → (O,S,P)；その鏡像異性体 (P,S,O) からは P が遊離
pick = [g for g in cand if g[0] == 'O']; assert pick == [('O', 'S', 'P')]
print('ORG-CBX-028 OK')
