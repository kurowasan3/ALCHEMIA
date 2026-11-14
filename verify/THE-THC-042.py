# 検証スクリプト THE-THC-042（グループ g12）  python3 -I この.py
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


it = check_item('/home/claude/chem6/outh/items/THE-THC-042/item.json',
                extra_eqs=['Zn + Cu^2+ -> Zn^2+ + Cu', 'Li + K+ -> Li+ + K'])
d = {'Zn': (130, 2639, -2046), 'Cu': (338, 2703, -2100), 'Li': (159, 520, -520), 'Na': (107, 496, -406), 'K': (89, 419, -322)}
def terms(a, b): return tuple(x - y for x, y in zip(d[a], d[b]))
t1 = terms('Zn', 'Cu'); print('Zn/Cu', t1, sum(t1)); assert t1 == (-208, -64, 54) and sum(t1) == -218
t2 = terms('Li', 'K'); print('Li/K', t2, sum(t2)); assert t2 == (70, 101, -198) and sum(t2) == -27
s = {m: sum(d[m]) for m in ('Li', 'Na', 'K')}; print(s); assert s == {'Li': 159, 'Na': 197, 'K': 186}
assert sorted(s, key=s.get) == ['Li', 'K', 'Na']
assert sorted(['Li', 'Na', 'K'], key=lambda m: d[m][1]) == ['K', 'Na', 'Li']
# 標準生成エンタルピー（NBS 表, kJ/mol）との照合：Zn2+(aq) -153.9, Cu2+(aq) +64.8, Li+ -278.5, Na+ -240.1, K+ -252.4
assert abs((-153.9 - 64.8) - sum(t1)) < 2
assert abs((-278.5 + 252.4) - sum(t2)) < 2
assert abs((s['Na'] - s['K']) - (-240.1 + 252.4)) < 2 and abs((s['K'] - s['Li']) - (-252.4 + 278.5)) < 2
# 標準電極電位の順 Li < K < Na（-3.04, -2.93, -2.71）と一致
E = {'Li': -3.04, 'K': -2.93, 'Na': -2.71}; assert sorted(E, key=E.get) == ['Li', 'K', 'Na']
print('ALL OK')
