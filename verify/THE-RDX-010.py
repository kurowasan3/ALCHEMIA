import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
E={'Cl2':1.36,'Br2':1.07,'Fe3+':0.77,'I2':0.54}
order=sorted(E, key=E.get, reverse=True); print(order); assert order==['Cl2','Br2','Fe3+','I2']
react=lambda ox, red_couple: E[ox]>E[red_couple]
assert react('Fe3+','I2') and not react('Fe3+','Br2') and react('Br2','Fe3+') and react('Cl2','Br2')
assert not react('I2','Fe3+') and react('Cl2','Fe3+') and react('Br2','I2')
from chemlib import check_equation
for eq in ['2Fe^2+ + Cl2 -> 2Fe^3+ + 2Cl-', 'Br2 + 2I- -> 2Br- + I2']: assert check_equation(eq) is None
assert not check_item('THE-RDX-010')
