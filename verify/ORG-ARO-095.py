# ORG-ARO-095 置換フェノールの酸の強さ（共鳴）と炭酸水素ナトリウムとの反応
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
mols = {'フェノール':'Oc1ccccc1','p-クレゾール':'Cc1ccc(O)cc1','m-ニトロフェノール':'Oc1cccc([N+](=O)[O-])c1',
        'p-ニトロフェノール':'Oc1ccc([N+](=O)[O-])cc1','2,4-ジニトロフェノール':'Oc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]',
        'ピクリン酸':'Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]'}
# 文献 pKa（25 ℃，水）: phenol 9.99, p-cresol 10.26, m-NP 8.36, p-NP 7.15, 2,4-DNP 4.09, picric 0.38; H2CO3 pKa1 6.35
pka = {'フェノール':10.0,'p-クレゾール':10.3,'m-ニトロフェノール':8.4,'p-ニトロフェノール':7.2,'2,4-ジニトロフェノール':4.1,'ピクリン酸':0.4}
for k,s in mols.items(): print(k, CalcMolFormula(Chem.MolFromSmiles(s)), pka[k], 'NaHCO3で CO2' if pka[k] < 6.4 else '')
for eq in ['C6H2(NO2)3OH + NaHCO3 -> C6H2(NO2)3ONa + H2O + CO2',
           'C6H3(NO2)2OH + NaHCO3 -> C6H3(NO2)2ONa + H2O + CO2',
           'C6H4(NO2)OH + NaOH -> C6H4(NO2)ONa + H2O']:
    print(eq, check_equation(eq))
# NaHCO3 との平衡定数の目安 K = 10^(6.4 - pKa)
for k in pka: print(k, 'K=10^', round(6.4-pka[k],1))
