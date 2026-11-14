# ORG-CBX-010 検証：C6H12O2 エステルの加水分解生成物の質量から酸・アルコール部分を決める
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
H, C, O = R(1), R(12), R(16)
pc = {'C': R('62.1'), 'H': R('10.3'), 'O': R('27.6')}
M = 116
n = {k: pc[k] / 100 * M / {'C': C, 'H': H, 'O': O}[k] for k in pc}
print({k: float(v) for k, v in n.items()}); assert [round(float(v)) for v in n.values()] == [6, 12, 2]
assert 6 * 12 + 12 + 2 * 16 == 116
nE = R(232, 100) / 116; assert nE == R(2, 100)
Mw = lambda g: R(g) / nE
assert (Mw('1.48'), Mw('1.20'), Mw('1.76'), Mw('0.92')) == (74, 60, 88, 46)
# 質量保存（水 0.0200 mol = 0.36 g）
assert R('2.32') + nE * 18 == R('1.48') + R('1.20') == R('1.76') + R('0.92')
def mw(s):
    m = Chem.MolFromSmiles(s); return sum({'C': 12, 'H': 1, 'O': 16}[a.GetSymbol()] * 1 + 0 for a in Chem.AddHs(m).GetAtoms())
# 分子量 46・60・74・88 の候補（高校範囲：飽和1価カルボン酸・飽和1価アルコール）
acids = {'HCOOH': 'OC=O', 'CH3COOH': 'CC(=O)O', 'C2H5COOH': 'CCC(=O)O', 'C3H7COOH': 'CCCC(=O)O'}
alcs = {'CH3OH': 'CO', 'C2H5OH': 'CCO', 'C3H7OH': 'CCCO', 'C4H9OH': 'CCCCO', 'C5H11OH': 'CCCCCO'}
pairs = [(a, b) for a, sa in acids.items() for b, sb in alcs.items() if mw(sa) + mw(sb) - 18 == 116]
print(pairs)
# P：問題文で「カルボン酸 1.48 g（M=74），アルコール 1.20 g（M=60）」と指定
P_pairs = [(a, b) for a, b in pairs if mw(acids[a]) == 74 and mw(alcs[b]) == 60]
assert P_pairs == [('C2H5COOH', 'C3H7OH')]
Q_pairs = [(a, b) for a, b in pairs if {mw(acids[a]), mw(alcs[b])} == {88, 46}]; print('Q', Q_pairs)
assert set(Q_pairs) == {('C3H7COOH', 'C2H5OH'), ('HCOOH', 'C5H11OH')}
# 構造
P, Q = 'CCC(=O)OC(C)C', 'CC(C)C(=O)OCC'
for s in (P, Q): assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C6H12O2'
iodo = Chem.MolFromSmarts('[CH3][CX4H1,CX4H2][OH]')
assert Chem.MolFromSmiles('CC(C)O').HasSubstructMatch(iodo) and not Chem.MolFromSmiles('CCCO').HasSubstructMatch(iodo)
for eq in ['C2H5COOCH(CH3)2 + H2O -> C2H5COOH + (CH3)2CHOH', '(CH3)2CHCOOC2H5 + H2O -> (CH3)2CHCOOH + C2H5OH']:
    assert check_equation(eq) is None, eq
print('ORG-CBX-010 OK')
