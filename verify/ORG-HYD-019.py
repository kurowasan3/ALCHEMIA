# ORG-HYD-019 検証：シン付加（冷・塩基性 KMnO4）/アンチ付加（Br2）の立体化学を 3D 座標で機械的に確かめる
import json, re, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def add_xy(smi, X='O', mode='syn'):
    """アルケン smi の C=C に X を2個付加。mode='syn' は同じ面（両面それぞれ），'anti' は反対面。生成物 SMILES の集合を返す"""
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=7); AllChem.MMFFOptimizeMolecule(m)
    b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); conf = m.GetConformer()
    P = lambda k: np.array(conf.GetAtomPosition(k))
    ni = [a.GetIdx() for a in m.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() != j]
    v1, v2 = P(j) - P(i), P(ni[0]) - P(i); n = np.cross(v1, v2); n /= np.linalg.norm(n)
    out = set()
    for s in (1, -1):
        rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.SINGLE)
        signs = (s, s) if mode == 'syn' else (s, -s)
        pos = []
        for k, sg in zip((i, j), signs):
            o = rw.AddAtom(Chem.Atom(8 if X == 'O' else 35)); rw.AddBond(k, o, Chem.BondType.SINGLE); pos.append((o, P(k) + sg * 1.4 * n))
        mm = rw.GetMol()
        c = mm.GetConformer()
        for o, p in pos: c.SetAtomPosition(o, p.tolist())
        for a in mm.GetAtoms(): a.SetNoImplicit(False)
        Chem.SanitizeMol(mm); Chem.AssignStereochemistryFrom3D(mm)
        out.add(Chem.MolToSmiles(Chem.RemoveHs(mm)))
    return out
mirror = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@', '!').replace('@', '@@').replace('!', '@')))
is_meso = lambda s: mirror(s) == Chem.MolToSmiles(Chem.MolFromSmiles(s)) and '@' in s
nst = lambda s: len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s))})

Z3, E3 = 'CC/C=C\\CC', 'CC/C=C/CC'
sZ, sE = add_xy(Z3), add_xy(E3)
print('Z-3-hexene syn', sZ, 'E-3-hexene syn', sE)
assert len(sZ) == 1 and is_meso(next(iter(sZ)))                     # (Z) → メソ体 1 種
assert len(sE) == 2 and all(not is_meso(x) for x in sE) and mirror(min(sE)) == max(sE)   # (E) → ラセミ体
assert nst('CCC(O)C(O)CC') == 3
bZ, bE = add_xy(Z3, 'Br', 'anti'), add_xy(E3, 'Br', 'anti')
print('Br2 anti Z', bZ, 'E', bE)
assert len(bE) == 1 and is_meso(next(iter(bE))) and len(bZ) == 2      # アンチ付加は逆
cp = add_xy('C1=CCCC1'); print('cyclopentene', cp)
assert len(cp) == 1 and is_meso(next(iter(cp))) and nst('OC1CCCC1O') == 3
mcp = add_xy('CC1=CCCC1'); print('1-methylcyclopentene', mcp)
assert len(mcp) == 2 and mirror(min(mcp)) == max(mcp) and nst('CC1(O)CCCC1O') == 4
pZ, pE = add_xy('C/C=C\\CC'), add_xy('C/C=C/CC'); print('2-pentene Z', pZ, 'E', pE)
assert len(pZ) == 2 and len(pE) == 2 and not (pZ & pE) and nst('CC(O)C(O)CC') == 4
# 熱した硫酸酸性 KMnO4：3-ヘキセン → プロピオン酸のみ（C6H12 で 1 種類のカルボン酸だけを与えるのは 3-ヘキセン）
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
alkenes = ['C=CCCCC','CC=CCCC','CCC=CCC','C=C(C)CCC','CC=C(C)CC','CC(C)=CCC','C=CC(C)CC','C=CCC(C)C','CC=CC(C)C','C=C(CC)CC','C=CC(C)(C)C','C=C(C)C(C)C','CC(C)=C(C)C']
assert len({Chem.MolToSmiles(Chem.MolFromSmiles(a)) for a in alkenes}) == 13 and all(CalcMolFormula(Chem.MolFromSmiles(a)) == 'C6H12' for a in alkenes)
def cleave(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s)); b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); m.RemoveBond(i, j)
    for k in (i, j):
        o = m.AddAtom(Chem.Atom(8)); m.AddBond(k, o, Chem.BondType.DOUBLE)
    mm = m.GetMol(); Chem.SanitizeMol(mm)
    fr = []
    for f in Chem.GetMolFrags(mm, asMols=True):
        s2 = Chem.MolToSmiles(f)
        fr.append(s2)
    return fr
one = [a for a in alkenes if len(set(cleave(a))) == 1]
print('single carbonyl fragment', one); assert one == ['CCC=CCC', 'CC(C)=C(C)C']   # 後者はアセトン（ケトン）でカルボン酸ではない
it = json.load(open('/home/claude/chem/out2/items/ORG-HYD-019/item.json'))
n = 0
for t in (it['problem'], it['answer'], it['solution']):
    for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', t):
        if eq.count('->') == 1 and 'boxed' not in eq and 'R' not in eq:
            r = check_equation(eq); assert r is None, (eq, r); n += 1; print('eq ok', eq)
assert n >= 1
print('OK')
