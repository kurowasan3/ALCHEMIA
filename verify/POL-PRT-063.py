# POL-PRT-063 検証：ピログルタミル基と C 末端アミドをもつトリペプチド誘導体
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
M={'H':1,'C':12,'N':14,'O':16}
def mw(s):
    m=Chem.AddHs(Chem.MolFromSmiles(s)); return sum(M[a.GetSymbol()] for a in m.GetAtoms())
S={
 'Glu':'NC(CCC(=O)O)C(=O)O','Phe':'NC(Cc1ccccc1)C(=O)O','Gly':'NCC(=O)O',
 'P (pGlu-Phe-Gly-NH2)':'O=C1CCC(N1)C(=O)NC(Cc1ccccc1)C(=O)NCC(N)=O',
 'cyclo(Gln-Phe-Gly)':'NC(=O)CCC1NC(=O)CNC(=O)C(Cc2ccccc2)NC1=O',
 'cyclo(Glu-Phe-Gly) gamma-amide variant':'NC(=O)C1CCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)N1',
 'glycinamide':'NCC(N)=O','pGlu-Phe':'O=C1CCC(N1)C(=O)NC(Cc1ccccc1)C(=O)O',
 'P ring-opened (Glu-Phe-Gly-NH2)':'NC(CCC(=O)O)C(=O)NC(Cc1ccccc1)C(=O)NCC(N)=O',
 'alpha-lactam alternative (gamma-linked)':'O=C1NC1CCC(=O)NC(Cc1ccccc1)C(=O)NCC(N)=O',
}
for k,v in S.items(): print(f'{k:40s}', CalcMolFormula(Chem.MolFromSmiles(v)), mw(v))
assert CalcMolFormula(Chem.MolFromSmiles(S['P (pGlu-Phe-Gly-NH2)']))=='C16H20N4O4' and mw(S['P (pGlu-Phe-Gly-NH2)'])==332
assert mw(S['glycinamide'])==74 and mw(S['pGlu-Phe'])==276 and mw(S['P ring-opened (Glu-Phe-Gly-NH2)'])==350
# 縮合の数：Glu+Phe+Gly+NH3 - 4H2O
tot=mw(S['Glu'])+mw(S['Phe'])+mw(S['Gly'])+17-4*18; print('sum check', tot); assert tot==332
# 環の数（アミノ酸単位を頂点・アミド結合を辺とするグラフ）：辺3，頂点3 → 独立な環 3-3+1=1
print('ring count in unit graph =', 3-3+1)
# P のラクタム環の員数
m=Chem.MolFromSmiles(S['P (pGlu-Phe-Gly-NH2)'])
print('ring sizes in P', [len(r) for r in m.GetRingInfo().AtomRings()])
# 遊離の NH2（第一級アミン，アミド以外）と COOH の有無
amine=Chem.MolFromSmarts('[NX3;H2;!$(NC=O)]'); acid=Chem.MolFromSmarts('C(=O)[OX2H1]')
for k in ['P (pGlu-Phe-Gly-NH2)','cyclo(Gln-Phe-Gly)','glycinamide','pGlu-Phe']:
    mm=Chem.MolFromSmiles(S[k]); print(k,'free amine',mm.HasSubstructMatch(amine),'COOH',mm.HasSubstructMatch(acid))
