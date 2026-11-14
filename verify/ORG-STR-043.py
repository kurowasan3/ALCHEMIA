# 検証スクリプト（グループ g12）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-STR-043/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
assert it['field'] in FIELDS, it['field']

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from sympy import Rational as R
from rdkit import Chem
from collections import Counter
smi = {'ア':'O=Cc1ccccc1','イ':'OC(=O)c1ccccc1','ウ':'OC(=O)c1ccccc1O','エ':'Cc1ccccc1','オ':'OCc1ccccc1',
       'カ':'Oc1ccccc1','キ':'c1ccccc1','ク':'C=Cc1ccccc1','ケ':'C#C'}
def comp(s):
    m = Chem.AddHs(Chem.MolFromSmiles(s)); return Counter(a.GetSymbol() for a in m.GetAtoms())
r73 = [k for k,s in smi.items() if R(comp(s)['C'], R(comp(s)['H'],2)) == R(7,3)]
r21 = [k for k,s in smi.items() if R(comp(s)['C'], R(comp(s)['H'],2)) == 2]
print(r73, r21); assert r73 == ['ア','イ','ウ'] and r21 == ['カ','キ','ク','ケ']
for k in r73:
    c = comp(smi[k]); o2 = R(2*c['C'] + R(c['H'],2) - c['O'], 2); print(k, o2)
    if k == 'イ': assert o2 == R('0.150')/R('0.0200')
    else: assert o2 != R('7.5')

print('ALL OK')
