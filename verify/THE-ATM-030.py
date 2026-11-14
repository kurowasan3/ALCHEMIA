# THE-ATM-030 verification: nuclear spin, unpaired electrons, ESR frequency
nuc = {'1H': (1, 1), '2H': (1, 2), '12C': (6, 12), '13C': (6, 13), '14N': (7, 14), '16O': (8, 16),
       '19F': (9, 19), '32S': (16, 32), '35Cl': (17, 35)}
nospin = [k for k, (z, a) in nuc.items() if z % 2 == 0 and a % 2 == 0]
print(nospin); assert nospin == ['12C', '16O', '32S']
Z = dict(H=1, C=6, N=7, O=8, Cl=17, Cu=29, Zn=30, Fe=26)
spec = {'CH3': Z['C'] + 3, 'NO': 15, 'NO2': 7 + 16, 'O3': 24, 'ClO2': 17 + 16, 'CN-': 6 + 7 + 1,
        'Cu2+': 27, 'Zn2+': 28, 'Fe3+': 23}
odd = [k for k, v in spec.items() if v % 2]
print(spec, odd); assert odd == ['CH3', 'NO', 'NO2', 'ClO2', 'Cu2+', 'Fe3+']
g, mu, h, B, c = 2.00, 9.27e-24, 6.63e-34, 0.330, 3.00e8
nu = g * mu * B / h; lam = c / nu
print(nu, lam, c / 9.23e9)
assert round(nu / 1e9, 2) == 9.23 and round(lam * 100, 2) == 3.25 and round(c / 9.23e9 * 100, 2) == 3.25
print('OK')
