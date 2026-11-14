ITEM = '/home/claude/chem5/out/items/POL-SUG-035/item.json'
# 検証スクリプト（グループ g09）  python3 -I この.py
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
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
assert [int(m.group(2)) for m in re.finditer(r'(^|\s)\((\d+)\)', it['answer'])] == list(range(1, np_+1))

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
for e in ['C6H10O5 + CS2 + NaOH -> C7H9O5S2Na + H2O','2C7H9O5S2Na + H2SO4 -> 2C6H10O5 + 2CS2 + Na2SO4','Cu(OH)2 + 4NH3 -> [Cu(NH3)4]^2+ + 2OH-']:
    assert check_equation(e) is None
A={'H':1.0,'C':12,'O':16,'S':32}
Mg = 6*12+10*1.0+5*16; assert Mg==162
n = 81.0/Mg; cs2 = 12+2*32
m_cs2 = n/2*cs2; print(m_cs2); assert sig(m_cs2,3)==19.0
assert sig(n*Mg,3)==81.0
Mta = 12*12+16*1.0+8*16; assert Mta==288   # C6H7O2(OCOCH3)3 = C12H16O8
print(n*Mta); assert sig(n*Mta,3)==144
print('ALL OK')
