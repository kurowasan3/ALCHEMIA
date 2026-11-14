# ORG-ARO-074 検証：ナフタレンのスルホン化・アルカリ融解
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
for e in ['C10H8 + H2SO4 -> C10H7SO3H + H2O', 'C10H7SO3H + NaHCO3 -> C10H7SO3Na + H2O + CO2',
          'C10H7SO3Na + 2NaOH -> C10H7ONa + Na2SO3 + H2O', 'C10H7ONa + CO2 + H2O -> C10H7OH + NaHCO3',
          '2C10H8O + 23O2 -> 20CO2 + 8H2O']:
    assert check_equation(e) is None, e
m = Chem.MolFromSmiles('Oc1ccc2ccccc2c1'); print(CalcMolFormula(m), MolWt(m))
n = Q('12.8') / 128; nB = n * Q('0.60')
mB = nB * 144; mCO2 = nB * 10 * 44
print(float(mB), float(mCO2))
assert round(float(mB), 1) == 8.6 and round(float(mCO2)) == 26
print('OK')
