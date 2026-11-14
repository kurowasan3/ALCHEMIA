# 検証スクリプト THE-PER-015  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-PER-015/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
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
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
import json
E = {e['sym']: e for e in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/elements.json', encoding='utf-8'))['elements']}
given = {'Na':23.0,'Mg':24.3,'Al':27.0,'Cl':35.5,'Ar':39.9,'K':39.1,'Co':58.9,'Ni':58.7,'Zn':65.4,'Ge':72.6,'In':114.8,'Te':127.6,'I':126.9}
for k,v in given.items():
    assert abs(E[k]['mass']-v) < 0.06, (k, E[k]['mass'])
# 周期表上の位置の確認（Ga の上下左右）
assert (E['Al']['group'],E['Al']['period']) == (13,3) and (E['In']['group'],E['In']['period']) == (13,5)
assert (E['Zn']['group'],E['Zn']['period']) == (12,4) and (E['Ge']['group'],E['Ge']['period']) == (14,4)
avg = (27.0+114.8+65.4+72.6)/4; print('avg', avg); assert round(avg) == 70
print('上下', (27.0+114.8)/2, '左右', (65.4+72.6)/2)
assert abs((27.0+114.8)/2-70.9) < 1e-9 and abs((65.4+72.6)/2-69.0) < 1e-9
pairs = {'a':('Na','Mg'),'b':('Cl','Ar'),'c':('Ar','K'),'d':('Co','Ni'),'e':('Te','I')}
for k,(x,y) in pairs.items(): assert E[y]['z'] == E[x]['z']+1
rev = [k for k,(x,y) in pairs.items() if given[x] > given[y]]
assert rev == ['c','d','e'], rev
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
print('Ga2O3 式量', round(2*69.7+48,1))

print('ALL OK')
