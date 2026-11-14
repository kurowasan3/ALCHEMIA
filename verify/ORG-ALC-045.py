# ORG-ALC-045 検証：カニッツァロ反応の重水素標識（反応式・酸化数・生成比）
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['2C6H5CHO + NaOH -> C6H5CH2OH + C6H5COONa', 'HCHO + C6H5CHO + NaOH -> HCOONa + C6H5CH2OH']:
    print(e, check_equation(e) or 'OK')
# 酸化数（カルボニル炭素）: 結合する O は −2 相当（二重結合で +2），H は −1，C–C は 0
print('アルデヒドの C', +2 - 1, ' アルコールの C', +1 - 2, ' カルボキシラートの C', +3)
# 仮説B：還元される分子（H 型/D 型）＋ 移る水素（供与分子が H 型/D 型），各 1/2 で独立
dist = {}
for red, don in itertools.product('HD', 'HD'):
    key = ''.join(sorted(red + don))
    dist[key] = dist.get(key, 0) + Fr(1, 4)
print('ベンジルアルコールの C 上の (H,D):', dist)   # HH 1/4, DH 1/2, DD 1/4
# 仮説A（溶媒 D2O から H/D）：C6H5CHO 単独 in D2O → C に D が 1 個
print('仮説A in D2O: C6H5CHDOD，仮説B: C6H5CH2OD')
# 交差反応：DCDO が酸化（DCOO-），C6H5CHO が還元され D を受け取る → C6H5CHDOH
print('交差: C6H5CHDOH + DCOO-')
