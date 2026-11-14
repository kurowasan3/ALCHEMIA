# 検証スクリプト THE-ION-098  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/outh/items/THE-ION-098/item.json'
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
import sympy as sp
C = 1.2e-22; KCu = 6.0e-36; KZn = 2.0e-22
H1 = 0.30 + 2 * 0.050
S1 = C / H1**2; cu = KCu / S1
print('(1)', H1, cu, 0.050 * S1); assert sig(H1) == 0.40 and sig(cu) == 8.0e-15 and 0.050 * S1 < KZn
x = sp.symbols('x', positive=True)
xa = float(sp.nsolve((0.050 - x) * C - KZn * (2 * x + 1e-5)**2, x, 0.04))
xb = (-0.60 + math.sqrt(0.36 + 0.48)) / 8
print('(2)', xa, xb, xa / 0.05, -math.log10(2 * xa + 1e-5))
assert abs(xa - xb) / xb < 1e-3 and sig(xa / 0.05 * 100) == 79 and round(-math.log10(2 * xa + 1e-5), 1) == 1.1
# 沈殿開始の確認
assert 0.05 * C / (1e-5)**2 > KZn
z3 = KZn / (C / (1e-5)**2); print('(3)', z3); assert sig(z3) == 1.7e-10
H4 = 2.7e-5 * 0.30 / 0.44; S4 = C / H4**2; z4 = KZn / S4
print('(4)', H4, -math.log10(H4), z4); assert round(-math.log10(H4), 1) == 4.7 and sig(z4) == 5.6e-10 and z4 < 1e-3 * 0.05
assert abs(math.log10(1.84) - 0.265) < 5e-4 and abs(math.log10(7.91) - 0.898) < 5e-4
assert abs(-math.log10(2.7e-5 * 0.20 / 0.54) - 5.0) < 0.01   # 緩衝液の初めの pH 5.0
for e in ['Zn^2+ + H2S -> ZnS + 2H+', 'Cu^2+ + H2S -> CuS + 2H+', 'CH3COO- + H+ -> CH3COOH', 'H2S <=> 2H+ + S^2-']:
    assert check_equation(e) is None, e
assert '79' in it['answer'] and '1.7\\times10^{-10}' in it['answer'] and '5.6\\times10^{-10}' in it['answer']

print('ALL OK')
