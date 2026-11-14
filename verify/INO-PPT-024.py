# INO-PPT-024 検証：試薬の組合せによる陽イオンの識別（ions.json から観察結果を生成して一意性を確認）
import sys, json, itertools
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
D = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))
R = D['results']; FLAME = {c['id']: c['flame'] for c in D['cations']}
def obs(ion, small, ex=None):
    a = R[ion][small]
    s = 'none' if a['kind'] in ('none', 'na') else a['color']
    if ex is None or s == 'none':
        return s
    return s + ('->dissolve' if R[ion][ex]['kind'] == 'dissolve' else '->remain')
cands = ['Ag+', 'Zn^2+', 'Al^3+', 'Pb^2+', 'Ba^2+', 'Ca^2+']
table = {i: (obs(i, 'NH3', 'NH3-ex'), obs(i, 'NaOH', 'NaOH-ex'), obs(i, 'Cl')) for i in cands}
for k, v in table.items(): print(k, v, FLAME[k])
given = {  # 問題文の表
 'A': ('白->dissolve', '白->dissolve', 'none'),
 'B': ('褐色->dissolve', '褐色->remain', '白'),
 'C': ('白->remain', '白->dissolve', '白'),
 'D': ('白->remain', '白->dissolve', 'none'),
 'E': ('none', 'none', 'none'),
}
ans = {'A': 'Zn^2+', 'B': 'Ag+', 'C': 'Pb^2+', 'D': 'Al^3+', 'E': 'Ba^2+'}
for L, row in given.items():
    m = [i for i in cands if table[i] == row]
    if L == 'E':
        m = [i for i in m if FLAME[i] == '黄緑']
    assert m == [ans[L]], (L, m)
# (3) アンモニア水だけで区別できない組
nh3 = {}
for i in ans.values(): nh3.setdefault(table[i][0], []).append(i)
pairs = [v for v in nh3.values() if len(v) > 1]
print('NH3 only, indistinct:', pairs); assert pairs == [['Pb^2+', 'Al^3+']]
# (4) Pb2+ と Al3+ を区別できる試薬
for reag in ['SO4', 'H2S-a', 'CrO4', 'NaOH-ex']:
    a, b = R['Pb^2+'][reag], R['Al^3+'][reag]
    print(reag, a.get('f'), a.get('color'), '|', b['kind'])
assert R['Pb^2+']['SO4']['kind'] == 'ppt' and R['Al^3+']['SO4']['kind'] == 'none'
assert R['Pb^2+']['H2S-a']['kind'] == 'ppt' and R['Al^3+']['H2S-a']['kind'] == 'none'
assert R['Pb^2+']['CrO4']['kind'] == 'ppt'   # Al3+ は na（高校教科書でクロム酸塩沈殿なし）
assert R['Pb^2+']['NaOH-ex']['kind'] == R['Al^3+']['NaOH-ex']['kind'] == 'dissolve'
assert check_equation('Ag2O + 4NH3 + H2O -> 2[Ag(NH3)2]+ + 2OH-') is None
print('INO-PPT-024 OK')
