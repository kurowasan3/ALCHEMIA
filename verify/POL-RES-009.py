# POL-RES-009 検証：試料の識別表
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from g12lib import check_item
S = {'ポリ塩化ビニル': dict(el={'C','H','Cl'}, sugar=False, arom=False),
     'ナイロン66': dict(el={'C','H','N','O'}, sugar=False, arom=False),
     '羊毛': dict(el={'C','H','N','O','S'}, sugar=False, arom=True),   # 芳香族アミノ酸を含むが繊維全体としてすすは多くない
     '綿': dict(el={'C','H','O'}, sugar=True, arom=False),
     'ポリスチレン': dict(el={'C','H'}, sugar=False, arom=True)}
def row(name):
    s = S[name]
    soot = name == 'ポリスチレン'   # 主鎖の繰り返し単位すべてにベンゼン環（C:H=1:1）
    return ('Cl' in s['el'], 'N' in s['el'], 'S' in s['el'], s['sugar'], soot)
table = {'A': (1,0,0,0,0), 'B': (0,1,0,0,0), 'C': (0,1,1,0,0), 'D': (0,0,0,1,0), 'E': (0,0,0,0,1)}
ans = {}
for k, t in table.items():
    m = [n for n in S if tuple(map(int, row(n))) == t]
    assert len(m) == 1, (k, m); ans[k] = m[0]
print(ans)
assert ans == {'A': 'ポリ塩化ビニル', 'B': 'ナイロン66', 'C': '羊毛', 'D': '綿', 'E': 'ポリスチレン'}
assert S['羊毛']['arom'] and not S['ナイロン66']['arom']   # キサントプロテインで区別
assert not check_item('POL-RES-009')
print('POL-RES-009 OK')
