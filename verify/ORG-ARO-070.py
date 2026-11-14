"""ORG-ARO-070 検証：反応式のつり合い・分子量・スズの質量"""
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

eqs = ['2C6H5NO2 + 3Sn + 14HCl -> 2C6H5NH3Cl + 3SnCl4 + 4H2O',
       'C6H5NH2 + NaNO2 + 2HCl -> C6H5N2Cl + NaCl + 2H2O',      # 15N を N として
       'C6H5N2+ + H2O -> C6H5OH + N2 + H+',
       'C6H5N2Cl + C6H5ONa -> C12H10N2O + NaCl']
for e in eqs:
    assert check_equation(e) is None, (e, check_equation(e))
azo = Chem.MolFromSmiles('Oc1ccc(cc1)N=Nc1ccccc1')
assert CalcMolFormula(azo) == 'C12H10N2O'
M = 12.00 * 12 + 1.00 * 10 + 14.00 + 15.00 + 16.00
assert abs(M - 199.00) < 1e-9
assert abs(14.00 + 15.00 - 29.00) < 1e-9
sn = 0.100 * 3 / 2 * 118.7
assert f'{sn:.3g}' == '17.8', sn
print('ORG-ARO-070 OK', M, sn)
