# THE-BND-114 verification: isomers of PFnCl(5-n) (trigonal bipyramid) and Berry pseudorotation paths
import itertools
from collections import deque
# positions 0,1 apical; 2,3,4 equatorial; D3h acting on positions
perms = set()
for k in range(3):
    for swap_ax in (False, True):
        for refl in (False, True):
            eq = [2 + (i + k) % 3 for i in range(3)]
            if refl: eq = [eq[0], eq[2], eq[1]]
            p = [1, 0] if swap_ax else [0, 1]
            perms.add(tuple(p + eq))
assert len(perms) == 12
def canon(t): return min(tuple(t[p[i]] for i in range(5)) for p in perms)
def isomers(nF):
    return {canon(tuple('F' if i in c else 'Cl' for i in range(5))) for c in itertools.combinations(range(5), nF)}
n = {k: len(isomers(k)) for k in range(6)}; print(n)
assert n[2] == 3 and n[1] == 2 and n[3] == 3
# kinds of F for each PF3Cl2 isomer (orbits of F positions under stabiliser)
for iso in sorted(isomers(3)):
    stab = [p for p in perms if tuple(iso[p[i]] for i in range(5)) == iso]
    orbits = {frozenset(p[i] for p in stab) for i in range(5) if iso[i] == 'F'}
    desc = sorted((('ax' if o <= {0, 1} else 'eq'), len(o)) for o in orbits)
    print(iso, desc)
    nax = (iso[0] == 'F') + (iso[1] == 'F')
    expect = {2: [('ax', 2), ('eq', 1)], 1: [('ax', 1), ('eq', 2)], 0: [('eq', 3)]}[nax]
    assert desc == expect
# pseudorotation on labelled ligands
def pr(st, piv):
    ax, eq = st
    return (tuple(sorted(x for x in eq if x != piv)), tuple(sorted([piv] + list(ax))))
def typ(st): return sum(x.startswith('F') for x in st[0])   # number of F apical
start = (('F1', 'F2'), ('Cl1', 'Cl2', 'Fs'))
assert typ(pr(start, 'Fs')) == 0 and typ(pr(start, 'Cl1')) == 1
def shortest(forbid0):
    seen = {start: []}; q = deque([start])
    while q:
        s = q.popleft()
        for p in s[1]:
            t = pr(s, p)
            if forbid0 and typ(t) == 0: continue
            if t not in seen:
                seen[t] = seen[s] + [p]; q.append(t)
    return min((len(v), v) for s, v in seen.items() if typ(s) == 2 and 'Fs' in s[0])
a, b = shortest(False), shortest(True); print(a, b)
assert a[0] == 2 and b[0] == 3
# check the explicit routes in the solution
s = pr(pr(start, 'Fs'), 'F1'); assert s == (('F2', 'Fs'), ('Cl1', 'Cl2', 'F1'))
s1 = pr(start, 'Cl1'); s2 = pr(s1, 'F1'); s3 = pr(s2, 'Cl2')
print(s1, s2, s3)
assert typ(s1) == 1 and typ(s2) == 1 and s3 == (('F1', 'Fs'), ('Cl1', 'Cl2', 'F2'))
print('OK')
