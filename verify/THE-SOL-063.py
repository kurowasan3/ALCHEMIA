# 検証スクリプト THE-SOL-063  python3 -I この.py
ITEM = '/home/claude/chem6/out/items/THE-SOL-063/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = set('原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split())
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
def check_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                r = check_equation(b)
                assert r is None, (b, r)
                print('式OK:', b); n += 1
    return n
def check_structure():
    assert it['field'] in FIELDS, it['field']
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
check_equations()
# ---- THE-SOL-063 の数値検証（SymPy・最小二乗）
import sympy as sp
c = [sp.Rational(x) for x in ('2.0', '4.0', '6.0', '8.0')]
Pi = [sp.Rational(x) for x in ('63.2', '132.8', '208.8', '291.2')]
r = [p/ci for p, ci in zip(Pi, c)]; print('Pi/c', [float(v) for v in r])
assert [float(v) for v in r] == [31.6, 33.2, 34.8, 36.4]
n = len(c); cx = sum(c)/n; ry = sum(r)/n
slope = sum((ci - cx)*(ri - ry) for ci, ri in zip(c, r))/sum((ci - cx)**2 for ci in c)
icpt = ry - slope*cx; print('slope', slope, 'intercept', icpt)
assert slope == sp.Rational(8, 10) and icpt == 30
assert all(ri == icpt + slope*ci for ci, ri in zip(c, r))   # 完全に直線
R = sp.Rational(83, 10)*10**3; T = 300
M = R*T/icpt; print('M', M); assert sig(float(M), 2) == 83000
M1 = R*T/r[-1]; print('M1', float(M1)); assert sig(float(M1), 2) == 68000 and M1 < M
PiY = icpt*8; print('Pi_Y', PiY); assert sig(float(PiY), 2) == 240
import xml.etree.ElementTree as ET
svg = open('/home/claude/chem6/out/items/THE-SOL-063/solution.svg', encoding='utf-8').read()
root = ET.fromstring(svg); assert root.get('viewBox')
assert not re.search(r'<script|foreignObject|\son\w+=|href=|base64', svg)
print('THE-SOL-063 OK')
