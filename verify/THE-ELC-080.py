# THE-ELC-080 検証：電気再生式脱イオン（EDI）の反応式とイオン移動
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
eqs=['H2O -> H+ + OH-','2H2O -> O2 + 4H+ + 4e-','2H2O + 2e- -> H2 + 2OH-','H+ + OH- -> H2O']
for e in eqs: print(e, check_equation(e))
# 樹脂の再生（一般式 R を含むので原子のつり合いを手で確認：R-SO3Na + H+ -> R-SO3H + Na+）
def cnt(side):
    from collections import Counter; return side
print('R-SO3Na + H+ -> R-SO3H + Na+ : 左 Na1 H1 電荷+1 / 右 H1 Na1 電荷+1')
print('R-N(CH3)3Cl + OH- -> R-N(CH3)3OH + Cl- : 左 Cl1 O1 H1 電荷-1 / 右 O1 H1 Cl1 電荷-1')
# イオンの移動：陽イオン→陰極(右)，陰イオン→陽極(左)
# 室II の Na+ は右へ → 膜Y は陽イオン交換膜, Cl- は左へ → 膜X は陰イオン交換膜
# 室I：Cl- が入り，陽極側から（電極反応で生じた H+ など）陽イオンは膜X を通れない → 濃縮側
print('膜X=陰イオン交換膜, 膜Y=陽イオン交換膜, 室II 脱塩')
