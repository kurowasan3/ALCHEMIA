ITEM = '/home/claude/chemx/x1/out/items/THE-ELC-079/item.json'
# 検証スクリプト（x1 g05）  python3 -I この.py
import sys, json, re, math
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS, it['field']
for f in it['format'].split(','):
    assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
import os
d = os.path.dirname(ITEM)
for k in ('problem_fig','solution_fig'):
    if it[k]: assert os.path.exists(os.path.join(d, it[k])), it[k]

# ---- 個別の検証 ----
from sympy import symbols, solve, Rational as R
a,b,c,dd = symbols('a b c d')
s = solve([R(1,2)*a-b, a-b-3*c, 2*a-2*b-4*c-2*dd], [b,c,dd])
assert [s[b].subs(a,6), s[c].subs(a,6), s[dd].subs(a,6)] == [3,1,1]
Li = {'L':R(1,2)*6, 'R':3}; Co={'L':6,'R':3+3}; O={'L':12,'R':6+4+2}
assert Li['L']==Li['R'] and Co['L']==Co['R'] and O['L']==O['R']
ox_before = 6*R(7,2); ox_after = 3*3 + 3*R(8,3); assert ox_before-ox_after == 4
M = 0.5*6.9+58.9+32; n = 18.9/M; chk('n', n, 0.200, 3)
chk('O2', n*0.5/6, 1.7e-2, 2)

print('ALL OK')
