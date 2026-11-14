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

it = check_item('/home/claude/chemx/kyodai/out/items/ORG-ARO-073/item.json')
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
import itertools
frags = ['C','CC','C=C','CCC','C=CC','CC=C','C(=C)C','C(C)C']   # 側鎖（先頭の原子が環に結合）
pos_sets = {1:[(0,)],2:[(0,1),(0,2),(0,3)],3:[(0,1,2),(0,1,3),(0,2,4)]}
found=set()
for k in (1,2,3):
    for fr in itertools.product(frags, repeat=k):
        for ps in pos_sets[k]:
            for perm in set(itertools.permutations(fr)):
                ring=['c']*6
                for p,f in zip(ps,perm): ring[p]='c('+f+')'
                smi=ring[0].replace('c','c1',1)+''.join(ring[1:5])+ring[5]+'1' if False else None
                atoms=[('c1' if i==0 else 'c') + (('('+perm[ps.index(i)]+')') if i in ps else '') for i in range(6)]
                smi=''.join(atoms)+''; smi=smi[:]  # 閉環
                smi=atoms[0]+''.join(atoms[1:5])+atoms[5].replace('c','c1',1) if False else ''.join(atoms[:5])+'c1'+(('('+perm[ps.index(5)]+')') if 5 in ps else '')
                m=Chem.MolFromSmiles(smi)
                if m is None: continue
                if D.CalcMolFormula(m)=='C9H10' and D.CalcNumRings(m)==1:
                    found.add(Chem.MolToSmiles(m))
print(found); assert len(found)==6
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
ct=[s for s in found if len(list(EnumerateStereoisomers(Chem.MolFromSmiles(s))))>1]
print('cis-trans:',ct); assert ct==['CC=Cc1ccccc1']
hyd=set()
for s in found:
    m=Chem.MolFromSmiles(s); rw=Chem.RWMol(m)
    for b in rw.GetBonds():
        if b.GetBondType()==Chem.BondType.DOUBLE and not b.GetIsAromatic(): b.SetBondType(Chem.BondType.SINGLE)
    hyd.add(Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rw))))
print(hyd); assert len(hyd)==5
import xml.dom.minidom; xml.dom.minidom.parse('/home/claude/chemx/kyodai/out/items/ORG-ARO-073/solution.svg')
print('OK')
