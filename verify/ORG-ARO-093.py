ITEM = '/home/claude/chemx/x2/outh/items/ORG-ARO-093/item.json'
# 検証スクリプト（x2 g03）  python3 -I この.py
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
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS
# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
# π電子数の数え方（規則どおり）と 4n+2 の判定
def is4n2(n): return n % 4 == 2
cnt = {'a': 2*2+2, 'b': 2*2+0, 'c': 3*2, 'd': 4*2, 'e': 2*2+2, 'f': 1*2}
assert cnt == {'a':6,'b':4,'c':6,'d':8,'e':6,'f':2}
arom = sorted(k for k,v in cnt.items() if is4n2(v)); assert arom == ['a','c','e','f']; print('芳香族', arom)
# RDKit の芳香族性の判定（ヒュッケル則）と照合
smi = {'a':'[cH-]1cccc1','b':'[CH+]1C=CC=C1','c':'[cH+]1cccccc1','d':'C1=CC=CC=CC=C1','e':'c1cc[nH]c1','f':'[cH+]1cc1'}
for k,s in smi.items():
    m = Chem.MolFromSmiles(s); ar = any(a.GetIsAromatic() for a in m.GetAtoms())
    assert ar == (k in arom), k; print(k, CalcMolFormula(m), 'RDKit芳香族' if ar else '非芳香族')
# トロポン（極性構造 6π）・シクロペンタジエノン（4π）・アズレン
for s, exp in [('O=C1C=CC=CC=C1', True), ('O=C1C=CC=C1', False), ('C1=CC=C2C=CC=C2C=C1', True)]:
    m = Chem.MolFromSmiles(s); ar = any(a.GetIsAromatic() for a in m.GetAtoms())
    assert ar == exp, s; print(s, CalcMolFormula(m), ar)
# アズレンの極性構造：五員環 C=C 2本+陰イオン=6，七員環 C=C 3本(陽イオン炭素0)=6；逆向き 4 と 8
assert is4n2(2*2+2) and is4n2(3*2) and not is4n2(2*2) and not is4n2(3*2+2)
assert is4n2(8+2)  # C8H8^2- 10π
print('ALL OK')
