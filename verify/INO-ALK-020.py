# INO-ALK-020 検証：カリウムの酸化物 P(K2O), Q(K2O2), R(KO2) の組成と反応
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from fractions import Fraction
K, O = 39.098, 15.999
def pct(nk, no): return 100*nk*K/(nk*K + no*O)
P, Q, R = pct(2,1), pct(2,2), pct(1,2)
print('K の質量%%: P %.2f  Q %.2f  R %.2f' % (P, Q, R))
assert round(P,1)==83.0 and round(Q,0)==71 and round(R,1)==55.0
# 問題文の値（3 桁）だけを使う（原子量を使わない）
p, q, r = 83.0, 71.0, 55.0
# (1) 単位格子：O は頂点 8×1/8 + 面心 6×1/2 = 4，K は内部 8 → K:O = 2:1
nO = Fraction(8,8) + Fraction(6,2); nK = 8
assert (nK, nO) == (8, 4)
# K 原子 1 個と O 原子 1 個の質量比
kr = (p/2)/(100-p)
for name, w in (('Q', q), ('R', r)):
    ratio = (w/(100-w))/kr      # nK/nO
    print(name, 'nK/nO =', round(ratio,4))
rq = (q/(100-q))/kr; rr = (r/(100-r))/kr
assert abs(rq-1) < 0.01 and abs(rr-0.5) < 0.01
# (3) 酸化数：O2^2- → -1，O2^- → -1/2
assert Fraction(-2,2) == -1 and Fraction(-1,2) == Fraction(-1,2)
# (4)(5) 反応式
for e in ['K2O2 + H2SO4 -> K2SO4 + H2O2', '4KO2 + 2CO2 -> 2K2CO3 + 3O2']:
    assert check_equation(e) is None, e; print('式OK', e)
print('(5) KO2 1 mol あたり CO2', Fraction(2,4), 'mol 吸収，O2', Fraction(3,4), 'mol 放出 → 気体は増える')
print('OK')
