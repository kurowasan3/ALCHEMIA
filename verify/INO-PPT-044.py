# 検証スクリプト INO-PPT-044  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/INO-PPT-044/item.json'
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
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
def rings(smi):
    m = Chem.MolFromSmiles(smi, sanitize=False)
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    ri = m.GetRingInfo()
    return sorted(len(r) for r in ri.AtomRings())
# 金属を [Co] とし，配位原子と単結合でつないだモデル
models = {
 'あ en': '[Co]1NCCN1',
 'い glycinate': '[Co]1NCC(=O)O1',
 'う oxalate': '[Co]1OC(=O)C(=O)O1',
 'え hexanediamine': '[Co]1NCCCCCCN1',
}
expect = {'あ en': [5], 'い glycinate': [5], 'う oxalate': [5], 'え hexanediamine': [9]}
for k, s in models.items():
    r = rings(s); print(k, r); assert r == expect[k]
# EDTA（6座）：SSSR の環の数と大きさ
edta = '[Co]12345OC(=O)CN1(CC(=O)O2)CCN3(CC(=O)O4)CC(=O)O5'
m = Chem.MolFromSmiles(edta, sanitize=False); m.UpdatePropertyCache(strict=False)
Chem.GetSymmSSSR(m)
r = sorted(len(x) for x in Chem.GetSymmSSSR(m))
print('EDTA', r); assert r == [5, 5, 5, 5, 5]
# (ii) en を展開してつり合いを確認
en = 'C2H8N2'
e2 = f'[Ni(H2O)6]^2+ + 3{en} -> [Ni({en})3]^2+ + 6H2O'
assert check_equation(e2) is None, check_equation(e2)
# 粒子数
assert (1 + 6, 1 + 6) == (7, 7) and (1 + 3, 1 + 6) == (4, 7)
assert '9 員環' in it['answer'] and '(3) b' in it['answer']

print('ALL OK')
