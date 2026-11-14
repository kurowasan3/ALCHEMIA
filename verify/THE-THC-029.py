"""THE-THC-029 検証：結合エネルギーによる孤立分子と網目（多面体）構造の安定性の比較（SymPy）"""
from sympy import Rational as R
E = {'C-O': 355, 'C=O': 805, 'Si-O': 450, 'Si=O': 630, 'N-N': 160, 'NtN': 945, 'P-P': 200, 'PtP': 490}
# 網目構造：XO2 1 単位あたり X-O 単結合 4 本（X は 4 本，O は 2 本で，各結合に X が 1 個）
# 反応「網目 XO2 -> 孤立分子 O=X=O」の ΔH = 切れる結合 − できる結合
dH_C = 4 * E['C-O'] - 2 * E['C=O']; dH_Si = 4 * E['Si-O'] - 2 * E['Si=O']
print('C:', dH_C, 'Si:', dH_Si); assert dH_C == -190 and dH_Si == 540
# 2X2 -> X4（正四面体：6 本の単結合）
dH_N = 2 * E['NtN'] - 6 * E['N-N']; dH_P = 2 * E['PtP'] - 6 * E['P-P']
print('N:', dH_N, 'P:', dH_P); assert dH_N == 930 and dH_P == -220
# 判定条件：多重度 k の結合 1 本が単結合 k 本より強いか
for a, b, k in (('C=O', 'C-O', 2), ('Si=O', 'Si-O', 2), ('NtN', 'N-N', 3), ('PtP', 'P-P', 3)):
    print(a, '/', b, '=', round(E[a] / E[b], 2), '比べる値', k, '多重結合が有利' if E[a] > k * E[b] else '単結合が有利')
assert E['C=O'] > 2 * E['C-O'] and E['Si=O'] < 2 * E['Si-O'] and E['NtN'] > 3 * E['N-N'] and E['PtP'] < 3 * E['P-P']
# 結合数の数え方の確認：網目 1 mol（X N_A 個）で X-O 結合 = 4 N_A（各 X が 4 本，各結合は X を 1 個だけ含む）
# P4：頂点 4，辺 6；N2 2 mol で三重結合 2 mol
assert 4 * 3 // 2 == 6
print('THE-THC-029 all OK')
