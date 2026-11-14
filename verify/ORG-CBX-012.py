"""ORG-CBX-012 検証：セッケン・合成洗剤の反応式のつり合い"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from g09lib import report_equations, formula
# 合成洗剤の陰イオン・生成物の分子式を RDKit でも確認
print('C12H25OSO3H', formula('CCCCCCCCCCCCOS(=O)(=O)O'))
print('C12H25C6H4SO3H', formula('CCCCCCCCCCCCc1ccc(S(=O)(=O)O)cc1'))
assert formula('CCCCCCCCCCCCOS(=O)(=O)O') == 'C12H26O4S'
assert formula('CCCCCCCCCCCCc1ccc(S(=O)(=O)O)cc1') == 'C18H30O3S'
report_equations('ORG-CBX-012')
