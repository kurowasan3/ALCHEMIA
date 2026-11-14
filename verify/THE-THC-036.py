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
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    print('小問数', np_, '答え', na)
    return it

it = check_item('/home/claude/chem6/out/items/THE-THC-036/item.json')
assert check_equation('CuSO4 + 5H2O -> CuSO4.5H2O') is None
m = 216
assert 16.0 + 200 == m and 25.0 + 191 == m and abs(5*0.100*18 - 9.0) < 1e-9
q1 = m*4.2*7.4; q2 = m*4.2*1.3
H1 = -q1/1000/0.100; H2 = q2/1000/0.100; Q = H1 - H2
print(H1, H2, Q); assert (round(H1), round(H2), round(Q)) == (-67, 12, -79)
# 吸湿の影響：16.0 g の試料のうち CuSO4 の物質量の割合 f が五水和物
H1t, H2t = H1, H2
for f in (0, 0.05, 0.1, 0.2):
    # 試料 16.0 g = n0(1-f)*160 + n0 f *250
    n0 = 16.0 / ((1-f)*160 + f*250)
    heat = -(n0*(1-f)*H1t + n0*f*H2t)   # 発熱量 kJ（正）
    H1m = -heat/0.100; Qm = H1m - H2t
    print(f, round(H1m,2), round(Qm,2))
    if f > 0: assert abs(Qm) < abs(Q)
print('THE-THC-036 OK')
