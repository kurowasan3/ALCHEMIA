# 検証スクリプト（グループ g06）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

FIELDS = set('原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split())

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
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
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    for f in it['format'].split(','):
        assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    assert it['answer'].count('$') % 2 == 0 and it['solution'].count('$') % 2 == 0 and it['problem'].count('$') % 2 == 0
    return it

def sig(x, n=2):
    from math import log10, floor
    if x == 0: return 0
    d = n - 1 - floor(log10(abs(x)))
    return round(x, d)

it = check_item('/home/claude/chem6/out/items/THE-ELC-046/item.json')
eta = {'Pt': 0.00, 'Cu': 0.35, 'Fe': 0.40, 'Zn': 1.00, 'Hg': 1.50}
def cathode(EH, Emetal, el):
    return 'H2' if EH - eta[el] > Emetal else 'metal'
assert cathode(0.00, -0.76, 'Pt') == 'H2' and cathode(0.00, -0.76, 'Zn') == 'metal'
assert cathode(-0.41, -1.80, 'Hg') == 'metal' and cathode(-0.41, -1.80, 'Fe') == 'H2'
assert abs((0.00-eta['Cu']) - (-0.76) - 0.41) < 1e-9 and abs(0.00 - (-0.76) - 0.76) < 1e-9
# 中性 (pH 7) の水素電極電位 -0.059*7 ≈ -0.41
assert abs(-0.0592*7 + 0.41) < 0.01
for e in ['Zn^2+ + 2e- -> Zn','2H+ + 2e- -> H2','2H2O + 2e- -> H2 + 2OH-','Na+ + e- -> Na','2Na + 2H2O -> 2NaOH + H2','Zn -> Zn^2+ + 2e-']:
    assert check_equation(e) is None
print('THE-ELC-046 OK')
