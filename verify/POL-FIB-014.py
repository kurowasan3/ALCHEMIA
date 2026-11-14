# 検証スクリプト POL-FIB-014  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/POL-FIB-014/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = "原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸".split()
assert it['field'] in FIELDS, it['field']
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
def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', sorted(eqs))
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    assert p.count('\\begin{enumerate}') <= 1
    n_items = len(re.findall(r'\\item\b', p)) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    assert n_items == len(nums), (n_items, nums)
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
from fractions import Fraction as F
def Xn(nA, nB, p, mono=F(0)):
    """nA: ジカルボン酸 mol, nB: ジアミン mol, p: 反応した COOH の割合（COOH が制限でない場合は NH2 で制限）"""
    before = nA + nB + mono
    cooh, nh2 = 2*nA + mono, 2*nB
    bonds = min(cooh*p, nh2)
    return before/(before - bonds)
x1 = Xn(F(1), F(1), F(995,1000)); x2 = Xn(F(1), F(1005,1000), F(1)); x3 = Xn(F(1), F(1005,1000), F(995,1000))
print(x1, x2, float(x3))
assert x1 == 200 and x2 == 401 and sig(float(x3), 3) == 134
r = 1/1.005; p = 0.995; car = (1+r)/(1+r-2*r*p); print('Carothers', car); assert abs(car - float(x3)) < 1e-9
# (4) 酢酸 q mol / 塩 1 mol（塩＝二酸1＋ジアミン1）：NH2 がすべて反応
q = F(2,199); x4 = Xn(F(1), F(1), F(1), mono=q); print(q, float(q), x4)
assert x4 == 200 and sig(float(q)*100, 2) == 1.0
assert Xn(F(1), F(1), F(1), mono=q*F(11,10)) < 200
print('ALL OK')
