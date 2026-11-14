# POL-PRT-056 検証：ペプチド Lys-Ala-Asp-Gly-Lys-Glu の電離する基
from rdkit import Chem
Q = Chem.MolFromSmiles('NCCCC[C@H](N)C(=O)N[C@@H](C)C(=O)N[C@@H](CC(=O)O)C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCC(=O)O)C(=O)O')
amine = Chem.MolFromSmarts('[NX3;H2][CX4]'); cooh = Chem.MolFromSmarts('C(=O)[OX2H1]')
na, nc = len(Q.GetSubstructMatches(amine)), len(Q.GetSubstructMatches(cooh))
print(na, nc); assert (na, nc) == (3, 3)  # 側鎖 2 + 末端 1 ずつ
assert (+na, -nc) == (3, -3) and na + nc == 6
assert na - nc == 0                 # pH 7
assert na - (nc - 1) == 1           # Asp→Asn
print('OK')
