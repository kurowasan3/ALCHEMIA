# POL-FIB-011 検証：(a) 確率の木を全探索して f(n) を厳密計算，(b) 漸化式と一致，(c) 長い鎖のモンテカルロで e^-2，(d) 質量計算・反応式
import sys, math, random
from fractions import Fraction
from functools import lru_cache
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

@lru_cache(None)
def exact(state):
    """state: タプル（1=未反応 OH）。期待される残りの OH 数を，全分岐を等確率でたどって求める"""
    pairs = [i for i in range(len(state)-1) if state[i] and state[i+1]]
    if not pairs:
        return Fraction(sum(state))
    tot = Fraction(0)
    for i in pairs:
        s = list(state); s[i] = s[i+1] = 0
        tot += exact(tuple(s))
    return tot / len(pairs)

f = {0: Fraction(0), 1: Fraction(1), 2: Fraction(0)}
for n in range(3, 12):
    f[n] = Fraction(2, n-1) * sum(f[k] for k in range(n-1))
for n in range(0, 10):
    assert exact(tuple([1]*n)) == f[n], n
print({n: str(f[n]) for n in range(7)})
assert f[3] == 1 and f[4] == Fraction(2, 3) and f[5] == 1 and f[6] == Fraction(16, 15)

# 長さ 4 で 0 個残る確率
@lru_cache(None)
def p_zero(state):
    pairs = [i for i in range(len(state)-1) if state[i] and state[i+1]]
    if not pairs:
        return Fraction(1 if sum(state) == 0 else 0)
    return sum(p_zero(tuple(0 if j in (i, i+1) else state[j] for j in range(len(state)))) for i in pairs) / len(pairs)
assert p_zero((1, 1, 1, 1)) == Fraction(2, 3)

# 長い鎖：漸化式で n=2000 の割合，およびモンテカルロ
fl = [0.0, 1.0, 0.0]; S = 0.0  # S = f(0)+...+f(n-2)
for n in range(3, 3001):
    S += fl[n-2]
    fl.append(2/(n-1)*S)
print('f(3000)/3000 =', fl[3000]/3000, 'e^-2 =', math.exp(-2))
assert abs(fl[3000]/3000 - math.exp(-2)) < 1e-3
random.seed(1)
N = 200000; free = [True]*N
order = list(range(N-1)); random.shuffle(order)
for i in order:   # 無作為な順に組を試す＝その時点の可能な組から一様に選ぶのと同等
    if free[i] and free[i+1]:
        free[i] = free[i+1] = False
print('MC', sum(free)/N); assert abs(sum(free)/N - 0.1353) < 0.003

# (3) 質量
m = 100 + 100/44*0.865/2*12
print('mass', m); assert f'{m:.3g}' == '112'
assert check_equation('CH3CH(OH)CH2CH(OH)CH3 + HCHO -> C6H12O2 + H2O') is None
assert CalcMolFormula(Chem.MolFromSmiles('CC1CC(C)OCO1')) == 'C6H12O2'
print('OK')
