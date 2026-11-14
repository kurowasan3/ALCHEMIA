# 検証スクリプト THE-SOL-062  python3 -I この.py
ITEM = '/home/claude/chem6/out/items/THE-SOL-062/item.json'
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
# ---- THE-SOL-062 の数値検証（SymPy）
import sympy as sp
Kf = sp.Rational(185, 100)
def mix(ions_a, ions_b, ppt):
    ions = {}
    for d in (ions_a, ions_b):
        for k, v in d.items(): ions[k] = ions.get(k, 0) + v
    for (c, an) in ppt:  # 1:1 沈殿・中和を可能な限り進める
        r = min(ions.get(c, 0), ions.get(an, 0))
        ions[c] -= r; ions[an] -= r
    return ions
q = sp.Rational(1, 100)
A = mix({'Ba': q, 'Cl': 2*q}, {'Na': 2*q, 'SO4': q}, [('Ba', 'SO4')])
B = mix({'Mg': q, 'Cl': 2*q}, {'Na': 2*q, 'SO4': q}, [])
C = mix({'Ba': q, 'OH': 2*q}, {'H': 2*q, 'SO4': q}, [('Ba', 'SO4'), ('H', 'OH')])
D = mix({'Ba': q, 'Cl': 2*q}, {'Na': 4*q, 'SO4': 2*q}, [('Ba', 'SO4')])
tot = {k: sum(v.values()) for k, v in zip('アイウエ', (A, B, C, D))}
print(tot)
assert tot == {'ア': sp.Rational(4,100), 'イ': sp.Rational(6,100), 'ウ': 0, 'エ': sp.Rational(7,100)}
fp = {k: -Kf*v/sp.Rational(2, 10) for k, v in tot.items()}
print({k: float(v) for k, v in fp.items()})
assert sorted(fp, key=lambda k: -fp[k]) == ['ウ', 'ア', 'イ', 'エ']
m = sp.Rational(925, 1000)/Kf; Wl = tot['ア']/m; ice = 200 - Wl*1000
print('ice', ice); assert sig(float(ice), 2) == 120
# (3) BaCl2 を s mol 加えたときの粒子数
def particles(s):
    ions = dict(D); ions['Ba'] = ions.get('Ba', 0) + s; ions['Cl'] += 2*s
    r = min(ions['Ba'], ions['SO4']); ions['Ba'] -= r; ions['SO4'] -= r
    return sum(ions.values())
h = sp.Rational(1, 1000)
s1 = (particles(5*h) - particles(0))/(5*h); s2 = (particles(20*h) - particles(15*h))/(5*h)
print('slopes', s1, s2); assert s1 == 1 and s2 == 3
assert particles(q) - particles(q - h) == h and particles(q + h) - particles(q) == 3*h
print('THE-SOL-062 OK')
