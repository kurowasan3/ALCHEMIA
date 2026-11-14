# 検証スクリプト POL-SUG-044  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/POL-SUG-044/item.json'
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
gal = 'OC[C@H]1O[C@@H](OC)[C@H](O)[C@@H](O)[C@H]1O'   # β-ガラクトシド（OC の C が結合先）
st = 'CCCCCCCCCCCCCCCCCC(=O)'
X = f'{st}OCC(OC(=O)CCCCCCCCCCCCCCCCC)CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O'
P = 'OCC(O)CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O'
R = 'O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O'
S = 'OCC(O)CO'
M = {'C': 12, 'H': 1, 'O': 16}
def mw(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    return sum(M[a.GetSymbol()] for a in m.GetAtoms())
for name, s, f, w in [('X', X, 'C45H86O10', 786), ('P', P, 'C9H18O8', 254), ('R', R, 'C6H12O6', 180), ('S', S, 'C3H8O3', 92)]:
    m = Chem.MolFromSmiles(s)
    print(name, CalcMolFormula(m), mw(s)); assert CalcMolFormula(m) == f and mw(s) == w
def nchiral(s):
    return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
print('不斉炭素 P', nchiral(P), 'S', nchiral(S)); assert nchiral(P) == 6 and nchiral(S) == 0
# 量的関係
n = 7.86 / 786
assert sig(2 * n * 40, 3) == 0.800 and sig(2 * n * 306, 3) == 6.12
assert mw('CCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]'.replace('.[Na+]', '')) + 23 == 306 + 1 - 1
assert '(3) R' in it['answer'] and 'P：6 個' in it['answer']

print('ALL OK')
