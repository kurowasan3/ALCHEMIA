# THE-BND-071 検証：等電子の組
v = {'H': 1, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7}
def info(atoms, q):
    return (len(atoms), sum(v[a] for a in atoms) - q)
ions = {'NO2+': (['N','O','O'], 1), 'N3-': (['N','N','N'], -1), 'NO2-': (['N','O','O'], -1),
        'NO3-': (['N','O','O','O'], -1), 'H3O+': (['O','H','H','H'], 1), 'NH2-': (['N','H','H'], -1)}
mols = {'CO2': ['C','O','O'], 'O3': ['O','O','O'], 'BF3': ['B','F','F','F'], 'NH3': ['N','H','H','H'], 'H2O': ['O','H','H'], 'CH4': ['C','H','H','H','H']}
exp = {'NO2+': 'CO2', 'N3-': 'CO2', 'NO2-': 'O3', 'NO3-': 'BF3', 'H3O+': 'NH3', 'NH2-': 'H2O'}
for k, (a, q) in ions.items():
    m = [n for n, b in mols.items() if info(b, 0) == info(a, q)]
    print(k, info(a, q), m); assert m == [exp[k]]
print('OK')
