# 検証スクリプト THE-GAS-057  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/outh/items/THE-GAS-057/item.json'
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
L = 44e3 / 18; c = 4.2; Lf = 6.0e3 / 18
m2 = 100 * c * 40 / L; print('(2)', m2); assert sig(m2) == 6.9
# 厳密：dm = m c dT / L
T = sp.symbols('T'); m = sp.Function('m')
sol = sp.dsolve(sp.Eq(m(T).diff(T), m(T) * c / L), m(T), ics={m(47): 100})
m_end = float(sol.rhs.subs(T, 7)); exact = 100 - m_end; print('(3) 厳密', exact); assert exact < m2
m47_0 = 100 * c * 47 / L; liq0 = 100 - m47_0
ice = liq0 * 44 / 50; print('(4)', m47_0, liq0, ice); assert sig(ice) == 81
# 0 ℃ の段階の数値計算（微小量ずつ蒸発させ，その熱を凝固でまかなう）
liq, ice_n, dy = liq0, 0.0, 1e-4
while liq > 0:
    dx = dy * L / Lf
    if dx + dy > liq:
        f = liq / (dx + dy); ice_n += dx * f; break
    liq -= dx + dy; ice_n += dx
print('数値', ice_n); assert abs(ice_n - ice) < 1e-2
assert '6.9' in it['answer'] and '81' in it['answer'] and '(3) 小さい' in it['answer']

print('ALL OK')
