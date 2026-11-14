# THE-BND-096 verification: Wiedemann-Franz law
data = {'Ag': (6.2e7, 429), 'Cu': (5.9e7, 401), 'Al': (3.7e7, 237)}
T = 300
L = {k: K / (s * T) for k, (s, K) in data.items()}
print(L)
assert [round(L[k] * 1e8, 1) for k in ('Ag', 'Cu', 'Al')] == [2.3, 2.3, 2.1]
Lm = sum(L.values()) / 3; KNa = Lm * 2.1e7 * T
print('Lmean', Lm, 'K_Na', KNa); assert round(KNa, -1) == 140
# T doubled, sigma halved -> K unchanged
assert abs(Lm * (2 * T) * (5.9e7 / 2) - Lm * T * 5.9e7) < 1e-9
print('OK')
