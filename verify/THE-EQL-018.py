import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item
# v1=k1[CO][H2]^2, v2=k2[CH3OH] を仮定し，平衡点から各操作直後の変化を調べる
co, h2, me = 1.0, 1.0, 1.0; k1 = 1.0; k2 = k1*co*h2**2/me
def rates(c, h, m, f1=1, f2=1): return k1*f1*c*h**2, k2*f2*m
v10, v20 = rates(co, h2, me)
def judge(v1, v2):
    s = lambda a, b: '増加' if a > b+1e-12 else ('減少' if a < b-1e-12 else '変化なし')
    d = '右' if v1 > v2+1e-12 else ('左' if v1 < v2-1e-12 else '移動しない')
    return s(v1, v10), s(v2, v20), d
res = [judge(*rates(2*co, 2*h2, 2*me)),
       judge(*rates(co, h2, me, f1=2.0, f2=3.0)),   # 逆反応(吸熱)の方が Ea 大→倍率大
       judge(*rates(co, h2, 0.5*me)),
       judge(*rates(co, h2, me)),
       judge(*rates(0.5*co, 0.5*h2, 0.5*me)),
       judge(*rates(co, h2, me, f1=5, f2=5))]
exp = [('増加','増加','右'),('増加','増加','左'),('変化なし','減少','右'),('変化なし','変化なし','移動しない'),('減少','減少','左'),('増加','増加','移動しない')]
for r, e in zip(res, exp): print(r, e); assert r == e
assert check_item('THE-EQL-018')
