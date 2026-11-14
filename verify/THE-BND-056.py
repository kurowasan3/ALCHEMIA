# 検証スクリプト THE-BND-056（グループ g12）  python3 -I この.py
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq): return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq): return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1: return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

FIELDS = set("""原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸""".split())

def check_item(path, extra_eqs=()):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            e = clean_ce(b)
            if not should_check(b): continue
            r = check_equation(e)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    for e in extra_eqs:
        r = check_equation(e); assert r is None, (e, r); n_eq += 1; print('式OK(追加):', e)
    print('検査した式の数:', n_eq)
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    for f in it['format'].split(','): assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it


it = check_item('/home/claude/chem6/out/items/THE-BND-056/item.json')
data = {'CCl4': (76.7, 29.8), 'C6H6': (80.1, 30.7), 'H2S': (-60.3, 18.7), 'CH3OCH3': (-24.8, 21.5), 'H2O': (100.0, 40.7), 'C2H5OH': (78.3, 38.6)}
r = {k: round(h*1000/(t+273.2)) for k, (t, h) in data.items()}; print(r)
assert r == {'CCl4': 85, 'C6H6': 87, 'H2S': 88, 'CH3OCH3': 87, 'H2O': 109, 'C2H5OH': 110}
normal = [r[k] for k in ('CCl4', 'C6H6', 'H2S', 'CH3OCH3')]
print('平均', sum(normal)/4); assert round(sum(normal)/4) == 87
h = 87*(-41+273.2)/1000; print('H2Se', h); assert round(h) == 20
T = 33200/87; t = T - 273.2; print('トルエン', T, t); assert round(t, -1) == 110
print('ALL OK')
