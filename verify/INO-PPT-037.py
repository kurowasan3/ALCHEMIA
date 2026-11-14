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

it = check_item('/home/claude/chem6/out/items/INO-PPT-037/item.json')
import itertools
cat = {'A':'Ba','B':'Ag','C':'Zn','D':'Cu'}   # (i)〜(iv) で決まる陽イオン
anions = ['Cl','NO3','SO4']
insol = {('Ag','Cl'):'AgCl', ('Ba','SO4'):'BaSO4'}   # 高校範囲で沈殿とする組（Ag2SO4 は図鑑 ions.json でも沈殿なし）
def ppts(x, y, an):
    s = set()
    for c, a in ((cat[x], an[y]), (cat[y], an[x])):
        if (c, a) in insol: s.add(insol[(c, a)])
    return s
def ok(an, use_vi=True):
    for k in an:
        if (cat[k], an[k]) in insol: return False          # 水溶性
    p = ppts('A','B',an)
    if p != {'AgCl'}: return False                         # (v) 白沈がすべて NH3 に溶ける
    if use_vi and ppts('A','C',an) != {'BaSO4'}: return False   # (vi)
    if ppts('A','D',an) or ppts('B','D',an): return False  # (vii)
    return True
sols = []; sols2 = []
for t in itertools.product(anions, repeat=4):
    an = dict(zip('ABCD', t))
    if ok(an): sols.append(t)
    if ok(an, False): sols2.append(t)
print('解:', sols); print('(vi) なし:', sols2)
assert sols == [('Cl','NO3','SO4','NO3')]
assert len(sols2) == 3
print('INO-PPT-037 OK')
