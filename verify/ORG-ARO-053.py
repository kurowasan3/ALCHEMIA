# 検証スクリプト ORG-ARO-053  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/ORG-ARO-053/item.json'
import sys, json, re
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
    assert it['field'] in FIELDS, it['field']
    assert all(f in ('qa','equation','compare','condition','pathway','calc','structure') for f in it['format'].split(','))
    import os
    for k in ('problem_fig','solution_fig'):
        if it[k]: assert os.path.exists(os.path.join(os.path.dirname(ITEM), it[k])), it[k]
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
import os
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
f = lambda s: CalcMolFormula(Chem.MolFromSmiles(s))
assert f('Nc1ccc(cc1)S(=O)(=O)O') == 'C6H7NO3S'
assert f('[NH3+]c1ccc(cc1)S(=O)(=O)[O-]') == 'C6H7NO3S'   # 双性イオンも同じ分子式
P = 'O=S(=O)([O-])c1ccc(cc1)/N=N/c1c(O)ccc2ccccc12.[Na+]'
Q = 'Oc1ccc2ccccc2c1/N=N/c1ccccc1'
print(f(P), f(Q))
assert f(P) == 'C16H11N2NaO4S' and f(Q) == 'C16H12N2O'
# カップリング位置：アゾ基の付いた炭素が OH の付いた炭素と隣接し，かつ環の縮合位の隣（1 位）
m = Chem.MolFromSmiles(Q)
azoC = [a for a in m.GetAtoms() if a.GetSymbol() == 'C' and any(n.GetSymbol() == 'N' for n in a.GetNeighbors()) and a.IsInRingSize(6) and m.GetRingInfo().NumAtomRings(a.GetIdx()) == 1 and any(n.GetSymbol() == 'C' and m.GetRingInfo().NumAtomRings(n.GetIdx()) == 2 for n in a.GetNeighbors())]
assert len(azoC) == 1
nbrs = azoC[0].GetNeighbors()
assert any(any(x.GetSymbol() == 'O' for x in n.GetNeighbors()) for n in nbrs)
print('カップリング位置（1 位）OK')
# 反応式の個別確認（HCl の内訳：塩に戻す 1 + ジアゾ化 2）
assert check_equation('NaO3SC6H4NH2 + HCl -> HO3SC6H4NH2 + NaCl') is None
assert check_equation('HO3SC6H4NH2 + NaNO2 + 2HCl -> HO3SC6H4N2Cl + NaCl + 2H2O') is None
for fig in ('problem.svg', 'solution.svg'):
    s = open(os.path.join(os.path.dirname(ITEM), fig), encoding='utf-8').read()
    assert 'viewBox' in s and 'script' not in s and 'foreignObject' not in s and 'currentColor' in s
print('図 OK')

print('ALL OK')
