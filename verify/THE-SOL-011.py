# 検証スクリプト（g13 作成）：python3 -I このファイル（作業ディレクトリのモジュールに依存しない）
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, Symbol, solve, Eq, nsimplify
ID = __file__.replace('\\', '/').split('/')[-1][:-3]
ITEM = json.load(open(f'/home/claude/chem/out/items/{ID}/item.json', encoding='utf-8'))
ARR = re.compile(r'<=>|->')


def ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0:
            return out
        k, d = j + 4, 1
        while d:
            d += {'{': 1, '}': -1}.get(s[k], 0)
            k += 1
        out.append(s[j + 4:k - 1])
        i = k


def check_all_equations():
    """問題文・答え・解説の \\ce{…} のうち矢印1つの完全な式をすべて検査"""
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(ITEM[key]):
            if len(ARR.findall(b)) != 1 or '\\boxed' in b or re.search(r'(?<![A-Za-z])[Rn](?![a-z])', b):
                continue
            r = check_equation(b)
            assert r is None, (key, b, r)
            n += 1
    print(f'反応式 {n} 本：原子数・電荷のつり合い OK')


def sig(x, n):
    """有効数字 n 桁に丸めた値"""
    return float(f'{float(x):.{n}g}')


def has(key, *frags):
    for f in frags:
        assert f in ITEM[key], (key, f)


def n_answers():
    a = re.findall(r'(?:^|(?<=[\s　]))\((\d+)\)', ITEM['answer'])
    k = ITEM['problem'].count('\\item')
    assert a == [str(i) for i in range(1, k + 1)], (a, k)
    print(f'小問 {k} 個と答えの番号が一致')

# ---------------------------------------------------------------- 本体
check_all_equations()
n_answers()
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors
smi = {'methanol': 'CO', '1-butanol': 'CCCCO', '1-hexanol': 'CCCCCCO', '1-octanol': 'CCCCCCCCO',
       'acetone': 'CC(C)=O', 'naphthalene': 'c1ccc2ccccc2c1', 'iodine': 'II',
       'sucrose': 'OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O'}
lp = {k: Crippen.MolLogP(Chem.MolFromSmiles(v)) for k, v in smi.items()}
for k, v in lp.items(): print(f'{k:12s} logP={v:+.2f}')
# 問題文の分子式とRDKitの分子式
for name, f in [('sucrose', 'C12H22O11'), ('acetone', 'C3H6O'), ('1-octanol', 'C8H18O'), ('naphthalene', 'C10H8')]:
    assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi[name])) == f, name
# (2) 水への溶けやすさ：logP の小さい順＝メタノール＞ブタノール＞ヘキサノール
assert lp['methanol'] < lp['1-butanol'] < lp['1-hexanol']
# (1) スクロースは親水的(logP<0)，1-オクタノール・ナフタレン・ヨウ素は疎水的(logP>1.5)，アセトンは中間
assert lp['sucrose'] < -2 and lp['1-octanol'] > 2 and lp['naphthalene'] > 2 and lp['iodine'] > 1
assert -1 < lp['acetone'] < 1
# 実測の水への溶解度（文献値，25 ℃，g/L）：メタノール 任意，1-ブタノール 約73，1-ヘキサノール 約5.9，1-オクタノール 約0.5
sol = {'1-butanol': 73, '1-hexanol': 5.9, '1-octanol': 0.5}
assert sol['1-butanol'] > sol['1-hexanol'] > sol['1-octanol']
# 図鑑：BaSO4 は白色沈殿（水に難溶）
ions = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json', encoding='utf-8'))
s = json.dumps(ions, ensure_ascii=False)
assert '"f": "BaSO4"' in s and '"kind": "ppt"' in s
# 三ヨウ化物イオンの式
assert check_equation('I2 + I- -> I3-') is None
has('answer', 'a：イ，b：ウ，c：ア，d：エ，e：イ', 'メタノール ＞ 1-ブタノール ＞ 1-ヘキサノール', '(4) 水')
print('THE-SOL-011 OK')
