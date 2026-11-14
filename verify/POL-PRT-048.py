# 検証スクリプト（グループ g08）  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/POL-PRT-048/item.json'
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if should_check(b):
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n_eq += 1
print('反応式のつり合い OK：', n_eq, '式')
# 小問数と答えの数の一致
m_en = re.search(r'\\begin\{enumerate\}(.*)\\end\{enumerate\}', it['problem'], re.S)
nq = m_en.group(1).count('\\item') if m_en else 0
assert '\\begin{itemize}' not in it['problem']
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問数', nq, '答えの数', na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


from rdkit import Chem
from itertools import combinations
side = {'グリシン':'[H]','アラニン':'C','セリン':'CO','システイン':'CS','リシン':'CCCCN','グルタミン酸':'CCC(=O)O','フェニルアラニン':'Cc1ccccc1','チロシン':'Cc1ccc(O)cc1'}
mw_aa = {'グリシン':75,'アラニン':89,'セリン':105,'システイン':121,'リシン':146,'グルタミン酸':147,'フェニルアラニン':165,'チロシン':181}
M = {'H':1,'C':12,'N':14,'O':16,'S':32}
def mw(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); return sum(M[a.GetSymbol()] for a in m.GetAtoms())
def aa(r):
    return 'NCC(=O)O' if r=='[H]' else f'NC({r})C(=O)O'
for k,r in side.items(): assert mw(aa(r)) == mw_aa[k], k
hits = []
for a,b in combinations(side, 2):
    ra = side[a] if side[a] != '[H]' else ''
    rb = side[b] if side[b] != '[H]' else ''
    smi = f"O=C1NC({ra})C(=O)NC1{rb}" if ra and rb else None
    if smi is None:
        r = ra or rb; smi = f"O=C1NC({r})C(=O)NC1" if r else 'O=C1NCC(=O)NC1'
    m = mw(smi)
    assert m == mw_aa[a] + mw_aa[b] - 36
    if m == 234: hits.append((a,b,smi))
print(hits)
assert sorted(h[:2] for h in hits) == sorted([('セリン','フェニルアラニン'),('アラニン','チロシン')])
arom = Chem.MolFromSmarts('c1ccccc1'); phen = Chem.MolFromSmarts('c[OX2H]'); S = Chem.MolFromSmarts('[#16]')
ok = []
for a,b,smi in hits:
    m = Chem.MolFromSmiles(smi)
    if m.HasSubstructMatch(arom) and not m.HasSubstructMatch(S) and not m.HasSubstructMatch(phen): ok.append((a,b))
print('条件を満たす組', ok); assert ok == [('セリン','フェニルアラニン')]

print('ALL OK')
