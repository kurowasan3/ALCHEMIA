# 検証スクリプト ORG-ALC-019  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/ORG-ALC-019/item.json'
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
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
    assert it['field'] in FIELDS, it['field']
    assert all(f in ('qa','equation','compare','condition','pathway','calc','structure') for f in it['format'].split(','))
    import os
    for k in ('problem_fig','solution_fig'):
        if it[k]: assert os.path.exists(os.path.join(os.path.dirname(ITEM), it[k])), it[k]
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
from sympy import Rational as R
H, C, O, Ca = R(10, 10), 12, 16, 40
acetone = 3*C + 6*H + O
ipa = 3*C + 8*H + O
caac = 4*C + 6*H + 4*O + Ca
assert acetone == 58 and ipa == 60 and caac == 158
ae_a = 3*acetone / (3*ipa + 294 + 4*98)
ae_b = acetone / caac
ae_c = acetone / ipa
benz, prop, o2, phenol = 6*C + 6*H, 3*C + 6*H, 2*O, 6*C + 6*H + O
ae_d1 = acetone / (benz + prop + o2)
ae_d2 = (phenol + acetone) / (benz + prop + o2)
vals = [sig(float(x*100), 3) for x in (ae_a, ae_b, ae_c, ae_d1)]
print('AE', vals, float(ae_d2*100))
assert vals == [20.1, 36.7, 96.7, 38.2] and ae_d2 == 1
yld = R(464, 10) / R(58)
frac = R(464, 10) / 60
print('yield', float(yld), 'mass frac', float(frac), 'ratio', float(frac/ae_a))
assert sig(float(yld*100), 3) == 80.0 and sig(float(frac*100), 3) == 77.3
assert abs(frac - ae_c*yld) < 1e-12
assert sig(float(frac/ae_a), 2) == 3.8
# 半反応式の確認
for e in ['CH3CH(OH)CH3 -> CH3COCH3 + 2H+ + 2e-', 'Cr2O7^2- + 14H+ + 6e- -> 2Cr^3+ + 7H2O']:
    assert check_equation(e) is None, e
# (4) 原子効率 100% の判定（生成物が 1 種類か）
eqs4 = {'ア': 'CH3CH=CH2 + H2O -> CH3CH(OH)CH3', 'イ': 'C2H5OH -> CH2=CH2 + H2O',
        'ウ': 'CH3COOH + C2H5OH -> CH3COOC2H5 + H2O', 'オ': 'C6H6 + HNO3 -> C6H5NO2 + H2O',
        'カ': '2CH2=CHCH3 + O2 -> 2CH3COCH3'}
full = [k for k, e in eqs4.items() if check_equation(e) is None and '+' not in e.split('->')[1]]
print('100%:', full)
assert full == ['ア', 'カ']

print('ALL OK')
