# 検証スクリプト POL-RUB-016  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/POL-RUB-016/item.json'
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
units = 5.40/54; n0 = 5.40/2.70e5; n1 = 5.40/5.40e4; O3 = 1.00e-4
print(units, n0, n1)
assert abs(O3/units - 1.0e-3) < 1e-15
cut = n1 - n0; f = cut/O3; hcho = O3 - cut
print(cut, f, hcho)
assert abs(f - 0.80) < 1e-9 and abs(hcho - 2.0e-5) < 1e-15
assert abs(n1/n0 - 5) < 1e-9
# 質量変化の見積り（主鎖切断で O2 原子分 +32，ビニル切断で +2 - HCHO 遊離分は高分子から失われる）
dm = cut*32 + hcho*(16-14)
print('mass change g', dm, dm/5.40); assert dm/5.40 < 1e-3
# (4) 反応性に差がある場合：主鎖 C=C の相対反応性 k
for ftrue in (0.5, 0.8):
    for k in (1, 2, 5):
        frac_main = ftrue*k/(ftrue*k + (1-ftrue))   # 少量のオゾンが主鎖に当たる割合
        print(ftrue, k, round(frac_main,4))
        if k > 1: assert frac_main > ftrue
        else: assert abs(frac_main-ftrue) < 1e-12
print('ALL OK')
