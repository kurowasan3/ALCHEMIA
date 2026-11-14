# POL-RUB-026 検証：N-イソプロピルアクリルアミドの合成式と分子式
# --- 共通：item.json の \ce{} 反応式のつり合い検査（単独で動く） ---
import json, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
def _ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k, d = j + 4, 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k
def check_item_equations(path):
    it = json.load(open(path, encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for ce in _ce_blocks(it[key]):
            arrows = re.findall(r'<=>|->|<-', ce)
            if len(arrows) != 1 or '\\boxed' in ce or '\\text' in ce or '\\mathrm' in ce:
                continue
            r = check_equation(ce)
            assert r is None, (key, ce, r)
            n += 1
    print('反応式', n, '本 OK')
    return it
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
check_item_equations('/home/claude/chemx/todai/out/items/POL-RUB-026/item.json')
A = Chem.MolFromSmiles('C=CC(=O)NC(C)C'); print('A', CalcMolFormula(A)); assert CalcMolFormula(A) == 'C6H11NO'
# 炭素数（疎水性の目安）：メチル 1 < エチル 2 < イソプロピル 3 → 白濁温度は逆順
order = sorted({'メチル体': 1, 'エチル体': 2, 'イソプロピル体': 3}.items(), key=lambda kv: -kv[1])
print([k for k, _ in order]); assert [k for k, _ in order] == ['イソプロピル体', 'エチル体', 'メチル体']
print('OK')
