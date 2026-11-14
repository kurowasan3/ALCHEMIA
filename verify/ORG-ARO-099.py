# ORG-ARO-099 閉環メタセシスと互変異性によるフェノールの合成：RDKit で生成物を機械的に求める
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
# 閉環メタセシス：分子内の2つの末端 =CH2 を除き（エチレンとして），内側の炭素どうしを C=C で結ぶ
rcm = AllChem.ReactionFromSmarts('([CH2:1]=[C:2].[CH2:3]=[C:4])>>([C:2]=[C:4]).[CH2:1]=[CH2:3]')
def RCM(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for ps in rcm.RunReactants((m,)):
        try:
            p = ps[0]; Chem.SanitizeMol(p); s = Chem.MolToSmiles(p)
            if p.GetRingInfo().NumRings() == 1: out.add(s)
        except Exception as e: pass
    return out
# ケト→芳香族フェノールへの互変異性：C(=O) の隣の sp3 CH から H を O へ（芳香環になる場合のみ）
taut = AllChem.ReactionFromSmarts('[C:1](=[O:2])[CX4;H1,H2:3]>>[c:1]([OH:2])[c:3]')
def to_phenol(smi):
    m = Chem.MolFromSmiles(smi)
    res=set()
    for ps in taut.RunReactants((m,)):
        p = ps[0]
        try:
            # 環の全原子を芳香族として再解釈：Kekule 化可能か
            s = Chem.MolToSmiles(p); q = Chem.MolFromSmiles(s)
            if q is not None and q.GetRingInfo().NumRings()==1 and all(a.GetIsAromatic() for a in q.GetAtoms() if a.IsInRing()):
                res.add(Chem.MolToSmiles(q))
        except Exception: pass
    return res
def phenol_from(diene):
    dien = RCM(diene); print('  diene', diene, CalcMolFormula(Chem.MolFromSmiles(diene)), '-> dienone', dien)
    ph=set()
    for d in dien:
        # 互変異性：dienone の異性化（H 移動）を、同じ分子式の芳香族フェノールとして判定
        m = Chem.MolFromSmiles(d)
        # 6員環で C=O 隣に sp3 CH があれば，その H を O に移すと芳香族
        for a in m.GetAtoms():
            if a.GetSymbol()=='C' and any(n.GetSymbol()=='O' and m.GetBondBetweenAtoms(a.GetIdx(),n.GetIdx()).GetBondTypeAsDouble()==2 for n in a.GetNeighbors()):
                ring_sp3 = [n for n in a.GetNeighbors() if n.IsInRing() and n.GetHybridization()==Chem.HybridizationType.SP3]
                print('   C=O neighbor sp3 ring C with H count:', [n.GetTotalNumHs() for n in ring_sp3])
        ph |= {can(x) for x in phenol_map.get(can(d), [])}
    return dien
# 期待されるフェノール（手で H を移したもの）と，dienone の対応を明示して確認
cases = {
 '(1) p-クレゾール': ('C=C(C)/C=C\\C(=O)CC=C', 'O=C1CC=C(C)C=C1', 'Cc1ccc(O)cc1'),
 '(2a) 2,5-ジメチルフェノール': ('C=C/C=C(/C)C(=O)CC(C)=C', 'O=C1CC(C)=CC=C1C', 'Cc1ccc(C)c(O)c1'),
 '(2b) 2,5-ジメチルフェノール': ('C=C/C(C)=C\\C(=O)C(C)C=C', 'O=C1C=C(C)C=CC1C', 'Cc1ccc(C)c(O)c1'),
 '(3) 6,6-ジメチル': ('C=C/C=C\\C(=O)C(C)(C)C=C', 'O=C1C=CC=CC1(C)C', None),
}
phenol_map = {}
for k,(diene, dienone, phenol) in cases.items():
    got = RCM(diene)
    ok = can(dienone) in got
    print(k, 'RCM ->', got, 'expected dienone match:', ok)
    if phenol:
        # dienone と phenol は同じ分子式（互変異性体）であること
        f1 = CalcMolFormula(Chem.MolFromSmiles(dienone)); f2 = CalcMolFormula(Chem.MolFromSmiles(phenol))
        print('   dienone', f1, 'phenol', f2, 'tautomer formula equal:', f1==f2)
        m = Chem.MolFromSmiles(dienone)
        sp3 = [a.GetTotalNumHs() for a in m.GetAtoms() if a.IsInRing() and a.GetHybridization()==Chem.HybridizationType.SP3]
        print('   ring sp3 carbons H counts:', sp3)
    else:
        m = Chem.MolFromSmiles(dienone)
        sp3 = [a.GetTotalNumHs() for a in m.GetAtoms() if a.IsInRing() and a.GetHybridization()==Chem.HybridizationType.SP3]
        print('   ring sp3 carbons H counts:', sp3, '-> no H: cannot aromatize')
print('2a and 2b give same phenol:', can(cases['(2a) 2,5-ジメチルフェノール'][2])==can(cases['(2b) 2,5-ジメチルフェノール'][2]))
# 2,5-ジメチルフェノールを与える前駆体を網羅：フェノールの2つの o 位のどちらを sp3 炭素(C6)とするか
# 置換基 R1..R5 を環の位置 2..6 として，C6 を o 位の一方に取る2通り
ph = Chem.MolFromSmiles('Cc1ccc(C)c(O)c1')
print('phenol formula', CalcMolFormula(ph))
# E/Z の確認：(1) の前駆体の幾何異性体
for smi in ('C=C(C)/C=C\\C(=O)CC=C', 'C=C(C)/C=C/C(=O)CC=C'):
    m = Chem.MolFromSmiles(smi)
    from rdkit.Chem import FindMolChiralCenters
    for b in m.GetBonds():
        if b.GetStereo()!=Chem.BondStereo.STEREONONE: print(smi, b.GetStereo())
# (4) 3-メチルフェノールの前駆体2種
for smi, exp in (('C=C/C(C)=C\\C(=O)CC=C','O=C1C=C(C)C=CC1'), ('C=C/C=C\\C(=O)CC(C)=C','O=C1CC(C)=CC=C1')):
    got = RCM(smi); print('(4)', smi, got, can(exp) in got, CalcMolFormula(Chem.MolFromSmiles(exp)), '-> 3-methylphenol', CalcMolFormula(Chem.MolFromSmiles('Cc1cccc(O)c1')))
