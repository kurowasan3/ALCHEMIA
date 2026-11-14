# 検証スクリプト THE-ION-074（グループ g12）  python3 -I この.py
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


it = check_item('/home/claude/chem6/out/items/THE-ION-074/item.json')
Ka = dict(AcOH=1.8e-5, H2CO3=4.5e-7, HCO3=4.7e-11, H3PO4=7.1e-3, H2PO4=6.3e-8, HPO4=4.5e-13, HCN=4.9e-10)
def K(a, b): return Ka[a]/Ka[b]
# (1) AcOH + PO4^3- -> HPO4^2-
k1 = K('AcOH', 'HPO4'); print('(1)', k1); assert k1 > 1e2
# (2) AcOH + HPO4^2- -> H2PO4-
k2 = K('AcOH', 'H2PO4'); print('(2)', k2); assert k2 > 1e2 and round(k2, -1) == 290
# (3) AcOH + H2PO4- -> H3PO4
k3 = K('AcOH', 'H3PO4'); print('(3)', k3); assert k3 < 1
# (4) H2CO3 + CN- ; その後 HCO3- + CN-
k4 = K('H2CO3', 'HCN'); k4b = K('HCO3', 'HCN'); print('(4)', k4, k4b); assert k4 > 1e2 and k4b < 1
# (5) H2PO4- + CO3^2- ; 続く反応はどれも起こらない
k5 = K('H2PO4', 'HCO3'); k5b = K('HPO4', 'HCO3'); k5c = K('H2PO4', 'H2CO3'); print('(5)', k5, k5b, k5c)
assert k5 > 1e2 and k5b < 1 and k5c < 1
# 等モルで混ぜるので，2段目以降は (1)(2) では酸が不足して起こらない（(1) の後に AcOH は残らない）
print('ALL OK')
