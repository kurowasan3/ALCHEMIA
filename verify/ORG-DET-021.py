"""ORG-DET-021 検証：溶媒の上下・抽出可否"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from g09lib import report_equations
d = {'ヘキサン': 0.655, '酢酸エチル': 0.902, 'クロロホルム': 1.48, 'アセトン': 0.784, 'メタノール': 0.792, '四塩化炭素': 1.59}
miscible = {'アセトン', 'メタノール'}
cls = {k: ('c' if k in miscible else ('a' if v < 1.0 else 'b')) for k, v in d.items()}
print(cls)
assert [k for k in d if cls[k] == 'a'] == ['ヘキサン', '酢酸エチル']
assert [k for k in d if cls[k] == 'b'] == ['クロロホルム', '四塩化炭素']
assert [k for k in d if cls[k] == 'c'] == ['アセトン', 'メタノール']
assert 1.33 > 1.0 > 0.713   # ジクロロメタンは下層（水層は上），エーテルは上層（水層は下）
report_equations('ORG-DET-021')
