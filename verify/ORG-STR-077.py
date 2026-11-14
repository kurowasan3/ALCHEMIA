# ORG-STR-077 検証：アルケン→クロロヒドリン（アンチ付加）→エポキシド（O⁻ と Cl がアンチの配座から背面攻撃）の立体化学
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def chlorohydrins():
    base='O[C{a}H](c1ccccc1)[C{b}H](Cl)C'
    out=[]
    for a in ('@','@@'):
        for b in ('@','@@'):
            s=Chem.MolToSmiles(Chem.MolFromSmiles(base.format(a=a,b=b)))
            out.append(s)
    return sorted(set(out))
ch=chlorohydrins(); print('クロロヒドリンの立体異性体',ch); assert len(ch)==4
def idx(m):
    # O, C1, C2, Cl, ipso, CH3
    patt=Chem.MolFromSmarts('[OX2H][CX4]([c:1])[CX4]([Cl])[CH3]')
    mt=m.GetSubstructMatch(Chem.MolFromSmarts('[OX2H]-[CX4](-c)-[CX4](-[Cl])-[CH3]'))
    O,C1,ip,C2,Cl,Me=mt
    return O,C1,C2,Cl,ip,Me
def epoxide_from(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m,randomSeed=7); AllChem.MMFFOptimizeMolecule(m)
    O,C1,C2,Cl,ip,Me=idx(m); conf=m.GetConformer()
    T.SetDihedralDeg(conf,O,C1,C2,Cl,180.0)
    d=T.GetDihedralDeg(conf,ip,C1,C2,Me)
    return 'cis' if abs(d)<90 else 'trans', d
def cip(smi):
    m=Chem.MolFromSmiles(smi)
    from rdkit.Chem import FindMolChiralCenters
    return FindMolChiralCenters(m)
res={}
for s in ch:
    g,d=epoxide_from(s); res[s]=g; print(s, cip(s), 'O と Cl をアンチにしたときの Ph–C–C–CH3 の二面角 %.0f'%d, '→', g)
assert sorted(res.values())==['cis','cis','trans','trans']
# アンチ付加：E 体・Z 体の 1-フェニルプロペンに，C2 の上側から Cl，C1 の下側から O が付く
def anti_add(alk):
    m=Chem.AddHs(Chem.MolFromSmiles(alk)); AllChem.EmbedMolecule(m,randomSeed=3); AllChem.MMFFOptimizeMolecule(m)
    mt=m.GetSubstructMatch(Chem.MolFromSmarts('c-[CH]=[CH]-[CH3]')); ip,C1,C2,Me=mt
    conf=m.GetConformer(); P=lambda i: np.array(conf.GetAtomPosition(i))
    nrm=np.cross(P(ip)-P(C1),P(C2)-P(C1)); nrm/=np.linalg.norm(nrm)
    em=Chem.RWMol(m)
    em.GetBondBetweenAtoms(C1,C2).SetBondType(Chem.BondType.SINGLE)
    o=em.AddAtom(Chem.Atom(8)); cl=em.AddAtom(Chem.Atom(17))
    em.AddBond(C1,o,Chem.BondType.SINGLE); em.AddBond(C2,cl,Chem.BondType.SINGLE)
    from rdkit.Geometry import Point3D
    c=em.GetConformer()
    po=P(C1)-1.4*nrm; pc=P(C2)+1.8*nrm
    c.SetAtomPosition(o,Point3D(*po)); c.SetAtomPosition(cl,Point3D(*pc))
    mm=em.GetMol()
    for a in mm.GetAtoms(): a.SetNoImplicit(False)
    Chem.SanitizeMol(mm)
    Chem.AssignStereochemistryFrom3D(mm)
    return Chem.MolToSmiles(Chem.RemoveHs(mm))
for alk,expect in (('C/C=C/c1ccccc1','trans'),('C/C=C\\c1ccccc1','cis')):
    p=anti_add(alk)
    can=Chem.MolToSmiles(Chem.MolFromSmiles(p))
    print(alk,'→ クロロヒドリン',can,cip(can),'→ エポキシド',res[can])
    assert res[can]==expect
# 鏡像のクロロヒドリンは鏡像のエポキシドを与えるので，4 種のクロロヒドリンは 4 種のエポキシド（cis 2 種・trans 2 種）に 1 対 1 で対応
ep=set()
for s in ['C[C@H]1O[C@@H]1c1ccccc1','C[C@@H]1O[C@H]1c1ccccc1','C[C@H]1O[C@H]1c1ccccc1','C[C@@H]1O[C@@H]1c1ccccc1']:
    ep.add(Chem.MolToSmiles(Chem.MolFromSmiles(s)))
print('エポキシドの立体異性体',ep); assert len(ep)==4
for s in ep: assert CalcMolFormula(Chem.MolFromSmiles(s))=='C9H10O'
print('OK')
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
assert check_equation('C6H5CH(OH)CH(Cl)CH3 + OH- -> C9H10O + Cl- + H2O') is None
assert check_equation('C6H5CH=CHCH3 + Cl2 + H2O -> C6H5CH(OH)CH(Cl)CH3 + HCl') is None
print('反応式 OK')
