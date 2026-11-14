# POL-PRT-064 検証：Gly・Ala・Phe 各2個の環状ヘキサペプチドの配列の数え上げと酵素切断による決定
from itertools import permutations
lin = set(permutations('GGAAFF'))
canon = lambda s: min(tuple(s[i:] + s[:i]) for i in range(len(s)))
cyc = sorted({canon(s) for s in lin})
print('鎖状', len(lin), '環状', len(cyc)); assert len(lin) == 90 and len(cyc) == 16
per3 = [c for c in cyc if c[:3] == c[3:]]; print('周期3', per3); assert len(per3) == 2
assert 6*(len(cyc) - len(per3)) + 3*len(per3) == 90
noadj = [c for c in cyc if all(c[i] != c[(i+1) % 6] for i in range(6))]
print('同じアミノ酸が隣り合わない', [''.join(c) for c in noadj]); assert len(noadj) == 5
def frags(c):   # Phe のカルボキシ基側（Phe-X 結合）だけを切る
    n = len(c); cut = [i for i in range(n) if c[i] == 'F']
    out = []
    for k, i in enumerate(cut):
        j = cut[(k+1) % len(cut)]
        f = []; p = (i+1) % n
        while True:
            f.append(c[p])
            if p == j: break
            p = (p+1) % n
        out.append(''.join(f))
    return sorted(out)
for c in noadj: print(''.join(c), frags(c))
P = [c for c in noadj if all(len(f) == 3 for f in frags(c)) and len(set(frags(c))) == 2]
print('P', P, frags(P[0])); assert len(P) == 1 and frags(P[0]) == ['AGF', 'GAF']
# P の 1 か所切断で得られる鎖状ペプチドの種類
nlin = lambda c: len({tuple(c[i:] + c[:i]) for i in range(6)})
print('P の1か所切断', nlin(P[0]), ' 周期3', [nlin(c) for c in per3]); assert nlin(P[0]) == 6 and all(nlin(c) == 3 for c in per3)
print('OK')
