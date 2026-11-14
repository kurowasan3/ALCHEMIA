# 検証スクリプト POL-FIB-018（グループ g10）  python3 -I この.py
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

FIELDS = set("原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸".split())

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
    return it

it = check_item('/home/claude/chem6/out/items/POL-FIB-018/item.json')
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def mw(smi):  # 原子量 H=1, C=12, O=16 で計算
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); w = {'H':1,'C':12,'O':16}
    return sum(w[a.GetSymbol()] for a in m.GetAtoms())
DMT = mw('COC(=O)c1ccc(cc1)C(=O)OC'); EG = mw('OCCO'); BHET = mw('OCCOC(=O)c1ccc(cc1)C(=O)OCCO')
unit = mw('OCCOC(=O)c1ccc(cc1)C=O') - 2  # -OCH2CH2O-CO-C6H4-CO- : C10H8O4
print(DMT, EG, BHET, unit); assert (DMT, EG, BHET, unit) == (194, 62, 254, 192)
from sympy import Rational as Q
n0 = Q(1940, 194); print('MeOH g', 2*n0*32); assert 2*n0*32 == 640
mol = n0 - Q(990,100); n = n0/mol; M = 192*n + 62
print('n', n, 'M', M, 'mass check', (n0*254 - Q(990,100)*62)/mol); assert n == 100 and M == 19262
assert (n0*254 - Q(990,100)*62)/mol == M
print('EG for n=200', n0 - n0/200); assert n0 - n0/200 == Q(995,100)
print('OK')
