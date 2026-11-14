CODE = 'ORG-STR-015'
# ---- 共通チェック（g10・難関編） ----
import sys, os, re, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
ITEM = os.path.join('/home/claude/chem/out2/items', CODE)
J = json.load(open(os.path.join(ITEM, 'item.json'), encoding='utf-8'))
FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
assert J['field'] in FIELDS, J['field']
assert J['difficulty'] in (4, 5) and 20 <= J['time_min'] <= 40
fm = J['format'].split(',')
assert 1 <= len(fm) <= 2 and all(f in 'qa equation compare condition pathway calc structure'.split() for f in fm)
assert 3 <= len(J['tags'].split(',')) <= 6
def _ce_all(s):
    out = []; i = 0
    while True:
        k = s.find('\\ce{', i)
        if k < 0: return out
        d = 0; j = k + 3
        while True:
            if s[j] == '{': d += 1
            elif s[j] == '}':
                d -= 1
                if d == 0: break
            j += 1
        out.append(s[k+4:j]); i = j
n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for e in _ce_all(J[key]):
        arrows = re.findall(r'<=>|->|<-', e)
        if len(arrows) != 1 or '\\boxed' in e or re.search(r'\bR\d?\b|\bn\b|\bAr\b|\bX\b', e) or '...' in e:
            continue
        r = check_equation(e)
        assert r is None, (key, e, r)
        n_eq += 1
n_items = J['problem'].count('\\item')
nums = re.findall(r'(?:^|\s)\((\d+)\)', J['answer'])
if n_items:
    assert [int(x) for x in nums] == list(range(1, n_items + 1)), (n_items, nums)
    assert J['problem'].count('\\begin{enumerate}') == 1
else:
    assert not nums, nums
for k in ('problem_fig', 'solution_fig'):
    if J[k]:
        sv = open(os.path.join(ITEM, J[k]), encoding='utf-8').read()
        assert 'viewBox' in sv and '<script' not in sv and 'foreignObject' not in sv and '<image' not in sv
        assert not re.search(r'\son[a-z]+=', sv) and 'href=' not in sv
for k in ('problem', 'answer', 'solution'):
    assert J[k].count('{') == J[k].count('}'), k
assert J['solution'].rstrip().endswith('\\hfill$\\square$')
print(f'[共通] {CODE}: 反応式 {n_eq} 本つり合いOK／小問 {n_items} 個と答えの番号一致／形式OK')
def formula(smi):
    return CalcMolFormula(Chem.MolFromSmiles(smi))
def mw(smi, table={'H':1.0,'C':12,'N':14,'O':16,'Cl':35.5,'Br':80,'S':32,'Na':23,'K':39}):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    return sum(table[a.GetSymbol()] for a in m.GetAtoms())
def canon(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))
# ---- ここから個別の検証 ----

import itertools
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
def rd(smi):
    out = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(smi), options=StereoEnumerationOptions(unique=True))}
    ach = 0
    for s in out:
        mir = s.replace('@@', '!').replace('@', '@@').replace('!', '@')
        ach += canon(mir) == s
    return len(out), ach
assert rd('ClC1CC=CCC1Cl') == (3, 1)
assert rd('ClC1CC(Cl)C(Cl)CC1Cl') == (7, 3)
assert rd('ClC1CCC(Cl)C(Cl)C1Cl') == (10, 2)
assert rd('ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl') == (9, 7)
# 独立な組合せ計算（平面六員環：回転・面内反転・裏返し）
OPS = []
for k in range(6):
    for refl in (False, True):
        for flip in (False, True):
            OPS.append(([((-i if refl else i) + k) % 6 for i in range(6)], flip, refl == flip))
def img(p, op):
    perm, flip, _ = op; q = [None]*6
    for i in range(6):
        if p[i] is not None:
            q[perm[i]] = ('D' if p[i] == 'U' else 'U') if flip else p[i]
    return tuple(q)
def classify(sub, cons=None):
    allp = set()
    for ud in itertools.product('UD', repeat=len(sub)):
        p = [None]*6
        for i, s in zip(sub, ud): p[i] = s
        allp.add(tuple(p))
    st, di, chi = set(), set(), {}
    for p in allp:
        if cons and not cons(p): continue
        K = lambda t: tuple(x or '-' for x in t)
        ks = min(K(img(p, o)) for o in OPS if o[2] and img(p, o) in allp); kd = min(K(img(p, o)) for o in OPS if img(p, o) in allp)
        st.add(ks); di.add(kd); chi[ks] = not any(img(p, o) == p for o in OPS if not o[2])
    return len(st), len(di), sum(not v for v in chi.values())
assert classify([0, 1, 3, 4]) == (7, 5, 3)
assert classify([0, 1, 2, 3]) == (10, 6, 2)
assert classify(list(range(6))) == (9, 8, 7)
assert classify([0, 1]) == (3, 2, 1)
assert classify([0, 1, 3, 4], lambda p: p[0] != p[1] and p[3] != p[4]) == (3, 2, 1)
assert classify([0, 1, 2, 3], lambda p: p[0] != p[1] and p[2] != p[3]) == (3, 2, 1)
def achiral(sub, pat):
    p = [None]*6
    for i, s in zip(sub, pat): p[i] = s
    p = tuple(p); return any(img(p, o) == p for o in OPS if not o[2])
def same(sub, a, b):
    pa = [None]*6; pb = [None]*6
    for i, s, t in zip(sub, a, b): pa[i] = s; pb[i] = t
    return any(img(tuple(pa), o) == tuple(pb) for o in OPS)
fig1245 = [('UUUU', True), ('UUDD', True), ('UUUD', False), ('UDDU', True), ('UDUD', False)]
fig1234 = [('UUUU', True), ('UUUD', False), ('UUDD', False), ('UUDU', False), ('UDDU', True), ('UDUD', False)]
for sub, figs in (([0, 1, 3, 4], fig1245), ([0, 1, 2, 3], fig1234)):
    for pat, ach in figs: assert achiral(sub, pat) == ach, (sub, pat)
    for a, b in itertools.combinations([p for p, _ in figs], 2): assert not same(sub, a, b), (a, b)   # 互いにジアステレオマー
assert same([0, 1, 3, 4], 'UUUD', 'UDUU')
print('ORG-STR-015 個別検証OK')
