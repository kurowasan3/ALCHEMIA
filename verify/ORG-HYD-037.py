# 検証スクリプト（グループ g07）  python3 -I この.py  （単独で動く）
import sys, json, re, os, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

FIELDS = set('''原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡
ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質
有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離
高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'''.split())
FORMATS = {'qa', 'equation', 'compare', 'condition', 'pathway', 'calc', 'structure'}

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

def check_item(path, dmin, dmax):
    it = json.load(open(path, encoding='utf-8'))
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    assert dmin <= it['difficulty'] <= dmax, it['difficulty']
    assert set(it['format'].split(',')) <= FORMATS, it['format']
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    d = os.path.dirname(path)
    for f in ('problem_fig', 'solution_fig'):
        if it[f]: assert os.path.exists(os.path.join(d, it[f])), it[f]
    assert '\\hfill$\\square$' in it['solution']
    return it

def sig(x, n):
    return float(f'{x:.{n}g}')
it = check_item('/home/claude/chemx/todai/out/items/ORG-HYD-037/item.json', 1, 3)
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
cases = {
 '1-methyl': ('CC1=CCCCC1', ['CC(=O)CCCCC(=O)O']),
 '3-methyl': ('CC1CCCC=C1', ['OC(=O)C(C)CCCC(=O)O']),
 '4-methyl': ('CC1CCC=CC1', ['OC(=O)CC(C)CCC(=O)O']),
 'methylene': ('C=C1CCCCC1', ['O=C1CCCCC1', 'O=C=O']),
}
def cnt(f, el):
    m = Chem.AddHs(Chem.MolFromSmiles(f)); return sum(1 for a in m.GetAtoms() if a.GetSymbol() == el)
res = {}
for k, (r, ps) in cases.items():
    dO = sum(cnt(p, 'O') for p in ps) - cnt(r, 'O'); dH = sum(cnt(p, 'H') for p in ps) - cnt(r, 'H')
    dC = sum(cnt(p, 'C') for p in ps) - cnt(r, 'C'); assert dC == 0
    e = 2*dO - dH
    chiral = any(Chem.FindMolChiralCenters(Chem.MolFromSmiles(p), includeUnassigned=True) for p in ps)
    res[k] = (CalcMolFormula(Chem.MolFromSmiles(r)), [CalcMolFormula(Chem.MolFromSmiles(p)) for p in ps], e, chiral)
    print(k, res[k])
assert res['1-methyl'][2] == 6 and all(res[k][2] == 8 for k in ('3-methyl', '4-methyl', 'methylene'))
assert [k for k in res if res[k][3]] == ['3-methyl', '4-methyl']
assert all(res[k][0] == 'C7H12' for k in res)
M = 7*12 + 12*1.0; n = 4.80/M; e = 6*n
avg = 0.2*3 + 0.8*5
mass = e/avg*(39.1 + 54.9 + 64.0); mass_all = e/5*158.0
print(M, n, e, avg, mass, mass_all)
assert sig(mass, 3) == 10.3 and sig(mass_all, 3) == 9.48
print('ORG-HYD-037 OK')
