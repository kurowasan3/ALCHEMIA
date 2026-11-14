# 検証スクリプト（京大 g04）  python3 -I この.py
import sys, json, re, math
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

FIELDS = ['原子の構造と電子配置','元素の周期表と周期性','化学結合と結晶','物質量と化学反応式','酸と塩基','酸化還元反応','電池と電気分解','気体の性質','溶液の性質','化学反応とエネルギー','反応速度と化学平衡','電離平衡と溶解平衡','ハロゲン・貴ガス・水素','酸素・硫黄','窒素・リン','炭素・ケイ素','アルカリ金属・アルカリ土類金属','アルミニウム・亜鉛・スズ・鉛','遷移元素（鉄・銅・銀・クロム・マンガン）','金属イオンの沈殿・錯イオン・系統分析','気体の製法と性質','有機化合物の構造と元素分析','脂肪族炭化水素','アルコール・アルデヒド・ケトン','カルボン酸・エステル・油脂','芳香族化合物','有機化合物の検出反応と分離','高分子の基礎と合成樹脂','合成繊維','ゴム・機能性高分子','糖類','アミノ酸・タンパク質・核酸']

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
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it

def sig(x, n=2):
    return float(f'{x:.{n-1}e}')

it = check_item('/home/claude/chemx/kyodai/out/items/ORG-CBX-052/item.json')
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
M = {'C':12,'H':1.0,'O':16}
def mw(f):
    c = {'C':0,'H':0,'O':0}
    for e,n in re.findall(r'([A-Z][a-z]?)(\d*)', f): c[e]+=int(n or 1)
    return sum(M[e]*n for e,n in c.items())
est = D.CalcMolFormula(Chem.MolFromSmiles('C=CCOC(=O)/C=C/c1ccccc1')); acid = D.CalcMolFormula(Chem.MolFromSmiles('OC(=O)C=Cc1ccccc1')); alc = D.CalcMolFormula(Chem.MolFromSmiles('C=CCO'))
assert (est,acid,alc)==('C12H12O2','C9H8O2','C3H6O'); assert mw(est)==188 and mw(acid)==148
n0 = 3.36/22.4/2; m0 = n0*188; w = 14.91-m0; nh = w/18
print(n0, m0, w, nh, nh/n0, nh*148, n0*5*22.4)
assert abs(n0-0.0750)<1e-12 and round(m0,1)==14.1 and round(nh/n0*100)==60 and sig(nh*148)==6.7 and sig(n0*5*22.4)==8.4
print('OK')
