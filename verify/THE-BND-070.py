# 検証スクリプト THE-BND-070  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/THE-BND-070/item.json'
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
import itertools
import sympy as sp
# 面心立方格子：格子点を数え，最近接の数を数える
pts = set()
for i, j, k in itertools.product(range(-2, 3), repeat=3):
    for d in [(0, 0, 0), (.5, .5, 0), (.5, 0, .5), (0, .5, .5)]:
        pts.add((i + d[0], j + d[1], k + d[2]))
dist = sorted(set(round(math.dist((0, 0, 0), p), 6) for p in pts if p != (0, 0, 0)))
dmin = dist[0]
n_nn = sum(1 for p in pts if p != (0, 0, 0) and abs(math.dist((0, 0, 0), p) - dmin) < 1e-6)
assert n_nn == 12 and abs(dmin - 1 / math.sqrt(2)) < 1e-5
# 単位格子中の数（頂点 1/8，面心 1/2）
assert 8 * sp.Rational(1, 8) + 6 * sp.Rational(1, 2) == 4
a = 0.562
d_nn = a / 1.41
print('最近接距離', d_nn, a / math.sqrt(2))
assert sig(d_nn) == 0.40 and sig(a / math.sqrt(2)) == 0.40
a2 = 1.41 * 0.276; print('a(仮)', a2); assert abs(a2 - 0.389) < 5e-4 and abs(0.389**3 - 0.0589) < 5e-5
rho = 4 * 18 / (6.0e23 * 0.0589e-21); print('仮の密度', rho, rho / 0.92)
assert sig(rho) == 2.0 and sig(rho / 0.92) == 2.2
assert '0.40' in it['answer'] and '2.0' in it['answer'] and '2.2 倍' in it['answer'] and '(4) い' in it['answer']

print('ALL OK')
