# 検証スクリプト POL-RUB-017  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/POL-RUB-017/item.json'
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
import random
random.seed(1)
def simulate(n0, nmol, p):
    N0 = n0*nmol; ncl = round(p*N0)
    while True:
        chosen = random.sample(range(N0), ncl)
        k = [0]*nmol
        for u in chosen: k[u//n0] += 1
        if min(k) >= 1: break
    eff = sum(x-1 for x in k); dang = 2*nmol
    return N0, eff, dang
for n0, nmol, p in [(100, 50, 0.05), (500, 20, 0.02)]:
    N0, eff, dang = simulate(n0, nmol, p)
    print(n0, p, eff, p*N0 - N0/n0, dang)
    assert eff == round(p*N0 - N0/n0)
f2 = (2/500)/(0.020 + 1/500); print('dangling', f2); assert sig(f2*100, 2) == 18
A = 0.010 - 1/2000; B = 0.012 - 1/250; print(A, B); assert A > B
p = A + 1/250; units = 100/(5*12+8*1.0); S = p*units/2*2*32
print(p, units, S); assert abs(p-0.0135) < 1e-12 and sig(S, 2) == 0.64
svg = open('/home/claude/chem3/outh/items/POL-RUB-017/solution.svg', encoding='utf-8').read()
assert 'viewBox' in svg and '<script' not in svg
print('ALL OK')
