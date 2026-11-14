# 検証スクリプト ORG-STR-033  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/ORG-STR-033/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
    assert it['field'] in FIELDS, it['field']
    assert all(f in ('qa','equation','compare','condition','pathway','calc','structure') for f in it['format'].split(','))
    import os
    for k in ('problem_fig','solution_fig'):
        if it[k]: assert os.path.exists(os.path.join(os.path.dirname(ITEM), it[k])), it[k]
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
import os
from itertools import product
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
# 図と同じジグザグ座標で molblock を作り，くさび（W）・破線（D）から立体を割り当てる
C = [(0, 0), (0.866, -0.5), (1.732, 0), (2.598, -0.5), (3.464, 0)]
O = {1: (0.866, -1.5), 2: (1.732, 1.0), 3: (2.598, -1.5)}
def mol_from(pattern, ohs=(1, 2, 3), ketone=()):
    atoms = [('C',) + c for c in C] + [('O',) + O[i] for i in ohs]
    bonds = [(1, 2, 1, 0), (2, 3, 1, 0), (3, 4, 1, 0), (4, 5, 1, 0)]
    for k, (i, w) in enumerate(zip(ohs, pattern)):
        if i in ketone:
            bonds.append((i + 1, 6 + k, 2, 0))
        else:
            bonds.append((i + 1, 6 + k, 1, 1 if w == 'W' else 6))
    L = ['', '  RDKit          2D', '', f'{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000']
    for a in atoms:
        L.append(f'{a[1]:10.4f}{a[2]:10.4f}{0:10.4f} {a[0]:<3} 0  0  0  0  0  0  0  0  0  0  0  0')
    for b in bonds:
        L.append(f'{b[0]:3d}{b[1]:3d}{b[2]:3d}{b[3]:3d}')
    L.append('M  END')
    m = Chem.MolFromMolBlock('\n'.join(L))
    return Chem.MolToSmiles(m)
def mirror(s):
    return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@', '!').replace('@', '@@').replace('!', '@')))
def groups(n, ohs):
    g = {}
    for pat in product('WD', repeat=n):
        g.setdefault(mol_from(pat, ohs), []).append(''.join(pat))
    return g
gd = groups(2, (1, 3))
print('diol', gd)
assert sorted(map(sorted, gd.values())) == [['DD', 'WW'], ['DW'], ['WD']]
meso_d = [v for s, v in gd.items() if mirror(s) == s]
assert meso_d == [['WW', 'DD']] or sorted(meso_d[0]) == ['DD', 'WW']
gt = groups(3, (1, 2, 3))
print('triol', gt)
assert sorted(map(sorted, gt.values())) == [['DDD', 'WWW'], ['DDW', 'DWW'], ['DWD', 'WDW'], ['WDD', 'WWD']]
key = {''.join(sorted(v)): s for s, v in gt.items()}
mesos = sorted(''.join(sorted(v)) for s, v in gt.items() if mirror(s) == s)
print('meso', mesos)
assert mesos == ['DDDWWW', 'DWDWDW']
assert mirror(key['WDDWWD']) == key['DDWDWW']   # 鏡像対
# EnumerateStereoisomers でも 4 種，擬似不斉 r/s はメソ体だけ
isos = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles('CC(O)C(O)C(O)C'), options=StereoEnumerationOptions(onlyUnassigned=False, unique=True))}
assert len(isos) == 4
for s in isos:
    cc = dict(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
    pseudo = any(v in ('r', 's') for v in cc.values())
    assert pseudo == (mirror(s) == s), (s, cc)
print('擬似不斉はメソ体のみ OK')
# (4) C3 の OH を除く
def deoxy(pat):
    return mol_from((pat[0], pat[2]), (1, 3))
assert deoxy('WWW') == deoxy('WDW') == mol_from('WW', (1, 3))
assert mirror(deoxy('WWW')) == deoxy('WWW')
assert deoxy('WWD') == deoxy('WDD') == mol_from('WD', (1, 3)) and mirror(deoxy('WWD')) != deoxy('WWD')
assert mirror(deoxy('WWD')) == deoxy('DWW')
# (5) 片端の酸化
def ox(pat, at):
    return mol_from(pat, (1, 2, 3), ketone=(at,))
a, b = ox('WWW', 1), ox('WWW', 3)
print('WWW', a, b)
assert a != b and mirror(a) == b
c, d = ox('WWD', 1), ox('WWD', 3)
print('WWD', c, d)
assert c != d and mirror(c) != d and mirror(c) != c and mirror(d) != d
s = open(os.path.join(os.path.dirname(ITEM), 'problem.svg'), encoding='utf-8').read()
assert 'viewBox' in s and 'script' not in s

print('ALL OK')
