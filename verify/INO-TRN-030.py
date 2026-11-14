# 検証スクリプト INO-TRN-030（グループ g07）  python3 -I この.py
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

def check_item(path, dmin, dmax):
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
    assert dmin <= it['difficulty'] <= dmax
    for f in it['format'].split(','):
        assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    nums = [int(m[1]) for m in re.findall(r'(^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, np_ + 1)), ('小問の数', np_, nums)
    assert '\\hfill$\\square$' in it['solution']
    return it

it = check_item('/home/claude/chem6/outh/items/INO-TRN-030/item.json', 4, 5)

H = {'FeS': -100, 'FeO': -272, 'Cu2S': -80, 'Cu2O': -170, 'SO2': -297}
ex = H['Cu2S'] + H['FeO'] - H['Cu2O'] - H['FeS']; assert ex == -82
fe = 2*H['FeO'] + 2*H['SO2'] - 2*H['FeS']; cu = 2*H['Cu2O'] + 2*H['SO2'] - 2*H['Cu2S']
print(fe, fe/3, cu, cu/3); assert fe == -938 and round(fe/3) == -313 and cu == -774 and cu/3 == -258
assert fe - cu == 2*ex
r3 = H['SO2'] - H['Cu2S'] - 2*H['Cu2O']; assert r3 == 123
tot = H['SO2'] - H['Cu2S']; assert tot == -217
from fractions import Fraction
x = Fraction(1, 1) / (1 + Fraction(1, 2)); assert x == Fraction(2, 3)
nCu = 1e6/63.5
o2, sio2, h2so4 = nCu*2.5*32/1e6, nCu*0.5*60/1e6, nCu*2*98/1e6
print(o2, sio2, h2so4); assert (round(o2,2), round(sio2,3), round(h2so4,2)) == (1.26, 0.472, 3.09)
assert check_equation('2CuFeS2 + 5O2 -> 2Cu + 2FeO + 4SO2') is None
assert check_equation('2FeO + SiO2 -> Fe2SiO4') is None

print('ALL OK')
