ITEM = '/home/claude/chemx/x1/out/items/ORG-ARO-089/item.json'
# 検証スクリプト（x1 g05）  python3 -I この.py
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS, it['field']
for f in it['format'].split(','):
    assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
import os
d = os.path.dirname(ITEM)
for k in ('problem_fig','solution_fig'):
    if it[k]: assert os.path.exists(os.path.join(d, it[k])), it[k]

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
def ncls(smi):
    m = Chem.MolFromSmiles(smi); r = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    return len({r[a.GetIdx()] for a in m.GetAtoms() if a.GetSymbol()=='C'})
res = [ncls(s) for s in ['OC(=O)c1ccccc1C(=O)O','OC(=O)c1cccc(C(=O)O)c1','OC(=O)c1ccc(C(=O)O)cc1']]
assert res == [4,5,3], res; print('classes', res)
def pos(p, k):  # ring positions 1..6
    return {((p-1+k) % 6)+1 for k in [k]}
def ortho(p): return {((p) % 6)+1, ((p-2) % 6)+1}
def para(p): return {((p+2) % 6)+1}
def meta(p): return {((p+1) % 6)+1, ((p-3) % 6)+1}
out = {}
for name, pm in [('o',2),('m',3),('p',4)]:
    occ = {1, pm}
    s = (meta(1)-occ) & ((ortho(pm)|para(pm))-occ)
    out[name] = s
print(out); assert out == {'o':{3,5}, 'm':set(), 'p':{3,5}}
assert F(Chem.MolFromSmiles('Cc1ccc(C(=O)O)cc1[N+](=O)[O-]')) == 'C8H7NO4'

print('ALL OK')
