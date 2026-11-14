# 検証スクリプト THE-EQL-080（グループ g07）  python3 -I この.py
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

def check_item(path, dmin, dmax):
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
    assert dmin <= it['difficulty'] <= dmax
    for f in it['format'].split(','):
        assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    nums = [int(m[1]) for m in re.findall(r'(^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, np_ + 1)), ('小問の数', np_, nums)
    assert '\\hfill$\\square$' in it['solution']
    return it

it = check_item('/home/claude/chem6/out/items/THE-EQL-080/item.json', 1, 3)

RT = 8.31e3*1000; V = 8.31
p2 = (0.100-0.020-0.050)*RT/V; print(p2); assert abs(p2-3.0e4) < 1e-6
assert 50 < p2 < 5.0e3 or (5.0e3 < p2 < 1e7)
ng = 5.0e3*V/RT; print(ng); assert abs(ng-0.0050) < 1e-12
caco3 = 0.100-0.020-ng; assert abs(caco3-0.075)<1e-12 and abs(0.090-caco3-0.015)<1e-12
assert 0.0 < 0.100-0.020-0.090 + 1 and 0.100-0.020-0.090 < 0   # CaO が余る

print('ALL OK')
