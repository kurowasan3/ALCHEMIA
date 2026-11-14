# 検証スクリプト ORG-ARO-047  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ARO-047/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b) and (' + ' in b or re.search(r'(^|\s)\d+\s*[A-Z(\[]', b)):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert len(nums) == n_items or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items, '答え数', len(nums))
check_structure()
check_all_equations()

# ---- 個別の検証 ----

from rdkit import Chem
# 環の位置 1..6 に置換基を置き，配向規則で選ばれる位置を求めて生成物を正規化 SMILES で数える
OP = {'OH': 3, 'CH3': 1}          # o,p 配向（値は活性化の強さ）
META = {'NO2', 'COOH'}
SM = {'OH': 'O', 'CH3': 'C', 'NO2': '[N+](=O)[O-]', 'COOH': 'C(=O)O', 'Br': 'Br', 'H': None}
def ring_smiles(subs):
    parts = []
    for i in range(1, 7):
        x = subs.get(i)
        parts.append('c' + ('1' if i == 1 else '') + (f'({SM[x]})' if x else ''))
    s = ''.join(parts)
    # 最後の原子で環を閉じる
    last = parts[-1]
    s = ''.join(parts[:-1]) + 'c1' + (last[1:] if len(last) > 1 else '')
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m)
def targets(subs):
    free = [i for i in range(1, 7) if i not in subs]
    want = {}
    for p, g in subs.items():
        if g in OP:
            pos = {(p % 6) + 1, ((p - 2) % 6) + 1, ((p + 2) % 6) + 1}
        elif g in META:
            pos = {((p + 1) % 6) + 1, ((p - 3) % 6) + 1}
        else:
            continue
        want[g] = {x for x in pos if x in free}
    return want
def main_positions(subs):
    want = targets(subs)
    sets = list(want.values())
    common = set.intersection(*sets)
    if common:
        return common, '一致'
    best = max((g for g in want if g in OP), key=lambda g: OP[g])
    return want[best], '競合'
def products(subs, new):
    pos, kind = main_positions(subs)
    out = set()
    for p in pos:
        d = dict(subs); d[p] = new
        out.add(ring_smiles(d))
    return out, kind
# (a) p-メチル安息香酸のニトロ化
pa, ka = products({1: 'CH3', 4: 'COOH'}, 'NO2')
assert ka == '一致' and pa == {Chem.MolToSmiles(Chem.MolFromSmiles('Cc1ccc(C(=O)O)cc1[N+](=O)[O-]'))}, pa
# (b) m-ニトロ安息香酸
pb, kb = products({1: 'COOH', 3: 'NO2'}, 'NO2')
assert kb == '一致' and pb == {Chem.MolToSmiles(Chem.MolFromSmiles('OC(=O)c1cc([N+](=O)[O-])cc([N+](=O)[O-])c1'))}, pb
# (c) p-クレゾールの臭素化
pc, kc = products({1: 'OH', 4: 'CH3'}, 'Br')
assert kc == '競合' and pc == {Chem.MolToSmiles(Chem.MolFromSmiles('Oc1ccc(C)cc1Br'))}, pc
# (3) ニトロトルエン 3 異性体
res = {}
for name, subs in {'o': {1: 'CH3', 2: 'NO2'}, 'm': {1: 'CH3', 3: 'NO2'}, 'p': {1: 'CH3', 4: 'NO2'}}.items():
    w = targets(subs)
    common = set.intersection(*w.values())
    prods = {ring_smiles({**subs, q: 'NO2'}) for q in common}
    res[name] = len(prods)
assert res == {'o': 2, 'm': 0, 'p': 1}, res
# (4) 2,4-DNT → TNT
w = targets({1: 'CH3', 2: 'NO2', 4: 'NO2'})
common = set.intersection(*w.values())
assert common == {6}, common
tnt = ring_smiles({1: 'CH3', 2: 'NO2', 4: 'NO2', 6: 'NO2'})
assert tnt == Chem.MolToSmiles(Chem.MolFromSmiles('Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]'))
print('配向の重ね合わせ OK', res)

print('ALL OK')
