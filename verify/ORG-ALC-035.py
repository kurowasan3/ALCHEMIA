# 検証スクリプト ORG-ALC-035  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/outh/items/ORG-ALC-035/item.json'
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
it = json.load(open(ITEM, encoding='utf-8'))


def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq


def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    return len(re.findall(r'<=>|->|<-', eq)) == 1


def check_common():
    assert it['field'] in FIELDS, it['field']
    for k in ('src_no', 'src_label', 'field', 'theme', 'difficulty', 'time_min', 'format', 'problem', 'answer',
              'solution', 'core', 'tags', 'verified', 'verify_note', 'problem_fig', 'solution_fig'):
        assert k in it, k
    for f in it['format'].split(','):
        assert f in ('qa', 'equation', 'compare', 'condition', 'pathway', 'calc', 'structure'), f
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == len(n_ans) and n_ans == [str(i + 1) for i in range(n_items)], (n_items, n_ans)
    eqs = []
    for key in ('problem', 'answer', 'solution'):
        for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[key]):
            if should_check(eq):
                e = check_equation(clean_ce(eq))
                assert e is None, (eq, e)
                eqs.append(eq)
    print('構造OK 小問', n_items, ' 検査した反応式', len(set(eqs)))
    for e in sorted(set(eqs)):
        print('   ', e)


def sig(x, n=2):
    if x == 0:
        return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))


check_common()

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def chir(m): return Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
cands = set()
for opos in [None, 0, 1, 2, 3]:
    atoms = ['O' if i == opos else 'C' for i in range(4)]
    for sub in ['C', 'O']:
        for sp in range(4):
            if atoms[sp] == 'O': continue
            a = [x if i != sp else x + '(' + sub + ')' for i, x in enumerate(atoms)]
            m = Chem.MolFromSmiles('c1ccc2c(c1)' + ''.join(a) + '2')
            if m and CalcMolFormula(m) == 'C10H12O':
                cands.add(Chem.MolToSmiles(m))
oh = Chem.MolFromSmarts('[OX2H]')
ethers = sorted(c for c in cands if not Chem.MolFromSmiles(c).HasSubstructMatch(oh))
print(len(cands), ethers); assert len(cands) == 8 and len(ethers) == 6
assert all(len(chir(Chem.MolFromSmiles(c))) == 1 for c in ethers)
names = {'2-メチルクロマン': 'CC1CCc2ccccc2O1', '3-メチルクロマン': 'CC1COc2ccccc2C1', '4-メチルクロマン': 'CC1CCOc2ccccc21',
         '1-メチルイソクロマン': 'CC1OCCc2ccccc21', '3-メチルイソクロマン': 'CC1Cc2ccccc2CO1', '4-メチルイソクロマン': 'CC1COCc2ccccc21'}
assert sorted(Chem.CanonSmiles(s) for s in names.values()) == ethers
hk = Chem.MolFromSmarts('[CX4]([OX2H])[OX2][#6]')
hemi = []
for k, s in names.items():
    m = Chem.MolFromSmiles(s)
    c = chir(m)[0][0]
    rw = Chem.RWMol(Chem.AddHs(m))
    h = [n.GetIdx() for n in rw.GetAtomWithIdx(c).GetNeighbors() if n.GetSymbol() == 'H'][0]
    rw.GetAtomWithIdx(h).SetAtomicNum(8)
    p = Chem.RemoveHs(rw.GetMol()); Chem.SanitizeMol(p)
    if p.HasSubstructMatch(hk): hemi.append(k)
    print(k, Chem.MolToSmiles(p), p.HasSubstructMatch(hk))
assert sorted(hemi) == sorted(['2-メチルクロマン', '1-メチルイソクロマン', '3-メチルイソクロマン'])
openf = {'2-メチルクロマン': 'CC(=O)CCc1ccccc1O', '1-メチルイソクロマン': 'CC(=O)c1ccccc1CCO', '3-メチルイソクロマン': 'CC(=O)Cc1ccccc1CO'}
phen = Chem.MolFromSmarts('c[OX2H]')
for k, s in openf.items():
    m = Chem.MolFromSmiles(s)
    assert CalcMolFormula(m) == 'C10H12O2' and not chir(m)
    print(k, s, 'phenol' if m.HasSubstructMatch(phen) else '')
assert [k for k, s in openf.items() if Chem.MolFromSmiles(s).HasSubstructMatch(phen)] == ['2-メチルクロマン']
# SVG の規格
import os
d = os.path.dirname(ITEM)
for f in ('problem.svg', 'solution.svg'):
    t = open(os.path.join(d, f)).read()
    assert 'viewBox' in t and '<script' not in t and 'foreignObject' not in t and 'href' not in t and not re.search(r'\son\w+=', t)

print('ALL OK')
