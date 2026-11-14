# 検証スクリプト（グループ g10）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/POL-FIB-016/item.json'
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
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS, it['field']
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
assert it['solution'].rstrip().endswith('\\hfill$\\square$')

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
M = {'H': 1.0, 'C': 12, 'N': 14, 'O': 16}
def fmass(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); tot = 0
    for a in m.GetAtoms():
        if a.GetSymbol() != '*': tot += M[a.GetSymbol()]
    return tot
# 繰り返し単位（* は鎖の続き）と アミド結合の数
units = {
 'ナイロン6':  ('*NCCCCCC(=O)*', 1, 'C6H11NO'),
 'ナイロン66': ('*NCCCCCCNC(=O)CCCCC(=O)*', 2, 'C12H22N2O2'),
 'ナイロン610':('*NCCCCCCNC(=O)CCCCCCCCC(=O)*', 2, 'C16H30N2O2'),
 'ナイロン11': ('*NCCCCCCCCCCC(=O)*', 1, 'C11H21NO'),
 'アラミド':   ('*Nc1ccc(cc1)NC(=O)c1ccc(cc1)C(=O)*', 2, 'C14H10N2O2'),
}
per = {}
for k, (smi, na, f) in units.items():
    m = Chem.MolFromSmiles(smi)
    fm = CalcMolFormula(m).replace('*', '')
    mass = fmass(smi)
    per[k] = mass / na
    print(k, CalcMolFormula(m), mass, per[k])
    assert CalcMolFormula(m).replace('*2', '') == f, (CalcMolFormula(m), f)
assert [round(per[k]) for k in units] == [113, 113, 141, 183, 119], per
v = {k: sig(100 / round(per[k]), 3) for k in ('ナイロン6', 'ナイロン610', 'ナイロン11')}
print(v)
assert v == {'ナイロン6': 0.885, 'ナイロン610': 0.709, 'ナイロン11': 0.546}
assert v['ナイロン6'] > v['ナイロン610'] > v['ナイロン11']
for s in ['0.885', '0.709', '0.546', '113 g', '141 g', '183 g', '119 g']:
    assert s in it['answer'], s
print('OK')
