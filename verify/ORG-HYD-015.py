"""ORG-HYD-015 検証：C5H10 アルケンの割り当ての一意性，ヨードホルム・銀鏡，シス-トランス，C6H12 の条件"""
import sys, itertools
sys.path.insert(0, '/home/claude/chem/work/g09')
from g09lib import report_equations, iodoform, aldehyde, has_cis_trans, can, formula
from alkenes import alkenes, hydrogenate, ozonolysis
HCHO = can('C=O')
al = alkenes(5)
print('C5H10 alkenes', al); assert len(al) == 5 and all(formula(a) == 'C5H10' for a in al)
info = {a: (hydrogenate(a), ozonolysis(a)) for a in al}
for a, v in info.items(): print(a, v)
def is_ketone(s): return not aldehyde(s)
sols = []
for P, Q, R, S in itertools.permutations(al, 4):
    T = info[P][0]
    if not (info[Q][0] == T and info[R][0] == T and info[S][0] != T): continue
    oP, oQ, oR, oS = (info[x][1] for x in (P, Q, R, S))
    if HCHO not in oP: continue
    K = [x for x in oP if x != HCHO]
    if len(K) != 1 or not is_ketone(K[0]): continue
    if HCHO not in oQ: continue
    L = [x for x in oQ if x != HCHO]
    if len(L) != 1 or not aldehyde(L[0]): continue
    if HCHO in oR or len(set(oR)) != 2: continue
    if sum(aldehyde(x) for x in oR) != 1: continue   # 実験3：M,N の一方だけ銀鏡
    if HCHO in oS or not all(aldehyde(x) for x in oS) or len(set(oS)) != 2: continue
    from rdkit import Chem
    nc = [Chem.MolFromSmiles(x).GetNumAtoms() for x in oS]
    if nc[0] == nc[1]: continue
    sols.append((P, Q, R, S, K[0], L[0], oR))
print('solutions', sols); assert len(sols) == 1
P, Q, R, S, K, L, oR = sols[0]
assert (P, Q, R, S) == (can('C=C(C)CC'), can('C=CC(C)C'), can('CC=C(C)C'), can('CC=CCC'))
N = [x for x in oR if aldehyde(x)][0]; M = [x for x in oR if not aldehyde(x)][0]
io = {k: iodoform(v) for k, v in dict(K=K, L=L, M=M, N=N).items()}
print('iodoform', io); assert io == dict(K=True, L=False, M=True, N=True)
ct = {k: has_cis_trans(v) for k, v in dict(P=P, Q=Q, R=R, S=S).items()}
print('cis-trans', ct); assert ct == dict(P=False, Q=False, R=False, S=True)
six = [a for a in alkenes(6) if len(set(ozonolysis(a))) == 1]
print('C6H12 single-product alkenes', six, [ozonolysis(a)[0] for a in six])
cond = [a for a in six if not aldehyde(ozonolysis(a)[0])]
print('no silver mirror', cond); assert cond == [can('CC(C)=C(C)C')]
report_equations('ORG-HYD-015')
