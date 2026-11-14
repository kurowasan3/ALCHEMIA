# THE-PER-019 検証：族の構成（第6周期まで）と金属元素の同定
# 高校の分類：非金属 = H, He, B, C, N, O, F, Ne, Si, P, S, Cl, Ar, As, Se, Br, Kr, Te, I, Xe, At, Rn（As, Te, At は非金属扱い）
nonmetal = set('H He B C N O F Ne Si P S Cl Ar As Se Br Kr Te I Xe At Rn'.split())
liquid = {'Hg', 'Br'}
gas = {'H', 'He', 'N', 'O', 'F', 'Ne', 'Cl', 'Ar', 'Kr', 'Xe', 'Rn'}
groups = {1: 'H Li Na K Rb Cs', 2: 'Be Mg Ca Sr Ba', 12: 'Zn Cd Hg', 13: 'B Al Ga In Tl',
          14: 'C Si Ge Sn Pb', 15: 'N P As Sb Bi', 16: 'O S Se Te Po', 17: 'F Cl Br I At', 18: 'He Ne Ar Kr Xe Rn'}
period = {}
for g, s in groups.items():
    for e in s.split():
        period[e] = {'H': 1, 'He': 1}.get(e, None)
order = {1: [1, 2, 3, 4, 5, 6], 2: [2, 3, 4, 5, 6], 12: [4, 5, 6], 13: [2, 3, 4, 5, 6], 14: [2, 3, 4, 5, 6], 15: [2, 3, 4, 5, 6], 16: [2, 3, 4, 5, 6], 17: [2, 3, 4, 5, 6], 18: [1, 2, 3, 4, 5, 6]}
for g, s in groups.items():
    for e, p in zip(s.split(), order[g]):
        period[e] = p
grp = {e: g for g, s in groups.items() for e in s.split()}
metals = [e for e in grp if e not in nonmetal]
def nm(g): return [e for e in groups[g].split() if e in nonmetal]
P = [e for e in metals if period[e] == 4 and len(nm(grp[e])) == 1 and nm(grp[e])[0] in gas]
Q = [e for e in metals if period[e] == 4 and any(x in liquid and x not in nonmetal for x in groups[grp[e]].split())]
R = [e for e in metals if period[e] == 3 and len(nm(grp[e])) == 1 and nm(grp[e])[0] not in gas and nm(grp[e])[0] not in liquid]
reacts_water = {'Li', 'Na', 'K', 'Rb', 'Cs', 'Ca', 'Sr', 'Ba'}
S = [e for e in metals if period[e] == 5 and len(nm(grp[e])) == 0 and e in reacts_water]
print(P, Q, R, S)
assert P == ['K'] and Q == ['Zn'] and R == ['Al'] and S == ['Sr']
# S 候補で水と反応しない Cd が除かれることを確認
print('S without water condition', [e for e in metals if period[e] == 5 and len(nm(grp[e])) == 0])
print('OK')
