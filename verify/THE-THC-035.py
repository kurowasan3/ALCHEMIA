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

it = check_item('/home/claude/chem6/out/items/THE-THC-035/item.json')
fH = {'H2O2(l)':-187.8,'H2O(l)':-285.8,'O2':0,'NH4Cl':-314.4,'NH3':-45.9,'HCl':-92.3,'N2':0,'H2':0,'O3':142.7}
S = {'H2O2(l)':109.6,'H2O(l)':69.9,'O2':205.2,'NH4Cl':94.6,'NH3':192.8,'HCl':186.9,'N2':191.6,'H2':130.7,'O3':238.9}
rx = {'ア':({'H2O2(l)':2},{'H2O(l)':2,'O2':1}), 'イ':({'NH4Cl':1},{'NH3':1,'HCl':1}), 'ウ':({'N2':1,'H2':3},{'NH3':2}), 'エ':({'O2':3},{'O3':2})}
expect = {'ア':'①','イ':'③','ウ':'②','エ':'④'}
for k,(L,R) in rx.items():
    dH = sum(fH[s]*c for s,c in R.items()) - sum(fH[s]*c for s,c in L.items())
    dS = sum(S[s]*c for s,c in R.items()) - sum(S[s]*c for s,c in L.items())
    cls = '①' if dH<0 and dS>0 else '②' if dH<0 and dS<0 else '③' if dH>0 and dS>0 else '④'
    print(k, round(dH,1), round(dS,1), cls); assert cls == expect[k]
print('THE-THC-035 OK')
