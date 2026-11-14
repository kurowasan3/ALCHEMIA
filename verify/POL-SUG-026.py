# 検証スクリプト POL-SUG-026  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/POL-SUG-026/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = "原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸".split()
assert it['field'] in FIELDS, it['field']
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
def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', sorted(eqs))
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    assert p.count('\\begin{enumerate}') <= 1
    n_items = len(re.findall(r'\\item\b', p)) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    assert n_items == len(nums), (n_items, nums)
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
SMI = {
 'beta-Glc': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O',
 'alpha-Glc': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O',
 'beta-Gal': 'C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O',
}
def stereo_diff(a, b):
    ca = Chem.FindMolChiralCenters(Chem.MolFromSmiles(a)); cb = Chem.FindMolChiralCenters(Chem.MolFromSmiles(b))
    return [(x, y) for x, y in zip(ca, cb) if x != y]
print('a/b', stereo_diff(SMI['beta-Glc'], SMI['alpha-Glc']))
print('Glc/Gal', stereo_diff(SMI['beta-Glc'], SMI['beta-Gal']))
assert len(stereo_diff(SMI['beta-Glc'], SMI['alpha-Glc'])) == 1
assert len(stereo_diff(SMI['beta-Glc'], SMI['beta-Gal'])) == 1
def axial_count(m, cid):
    conf = m.GetConformer(cid)
    ring = [r for r in m.GetRingInfo().AtomRings() if len(r) == 6][0]
    P = np.array([list(conf.GetAtomPosition(i)) for i in ring]); c = P.mean(0)
    nrm = np.linalg.svd(P - c)[2][2]
    cnt = 0
    for i in ring:
        a = m.GetAtomWithIdx(i)
        if a.GetSymbol() != 'C': continue
        for nb in a.GetNeighbors():
            if nb.GetIdx() in ring or nb.GetSymbol() == 'H': continue
            v = np.array(conf.GetAtomPosition(nb.GetIdx())) - np.array(conf.GetAtomPosition(i))
            cos = abs(v.dot(nrm))/np.linalg.norm(v)
            if cos > 0.75: cnt += 1
    return cnt
res = {}
for name, smi in SMI.items():
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    cids = AllChem.EmbedMultipleConfs(m, numConfs=300, randomSeed=7)
    out = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    E = [e for conv, e in out]
    data = sorted((E[k], axial_count(m, c)) for k, c in enumerate(cids))
    best = data[0]; counts = sorted(set(a for e, a in data))
    # 各 axial 数ごとの最低エネルギー
    emin = {}
    for e, a in data: emin[a] = min(emin.get(a, 1e9), e)
    res[name] = (best, emin); print(name, 'best', best, 'emin by axial', {k: round(v-best[0],1) for k,v in emin.items()})
assert res['beta-Glc'][0][1] == 0 and res['alpha-Glc'][0][1] == 1 and res['beta-Gal'][0][1] == 1
assert 5 in res['beta-Glc'][1] and res['beta-Glc'][1][5] > res['beta-Glc'][1][0]
print('ALL OK')
