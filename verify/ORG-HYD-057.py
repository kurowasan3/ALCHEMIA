ITEM = '/home/claude/chemk/out/items/ORG-HYD-057/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
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
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# --- 個別の検証 ---
h, c, NA = 6.63e-34, 3.00e8, 6.02e23
E = NA*h*c/475e-9/1000
print('E =', E); assert sig(E, 3) == 252
lam = {4:300,5:335,6:365,7:395,8:420,9:440,10:460,11:475}
def color(l):
    if l < 400: return '無色'
    for lo, hi, col in [(400,435,'黄緑'),(435,480,'黄'),(480,500,'橙〜赤'),(500,560,'赤紫'),(560,610,'青'),(610,750,'青緑')]:
        if lo <= l < hi: return col
assert color(lam[11]) == '黄'
# (2) 端に付加
assert lam[10] == 460 and color(lam[10]) == '黄'
# (3) k 番目に付加 → 最長の共役部分
ok = [k for k in range(1, 12) if max(k-1, 11-k) < 4 or lam.get(max(k-1, 11-k), 0) < 400]
print('無色になる位置', ok); assert ok == [4,5,6,7,8]
assert color(lam[8]) == '黄緑'
# (4)
assert color(lam[10]) == '黄' and color(lam[5]) == '無色'
# 表の単調性（共役が長いほど長波長）
v = [lam[k] for k in sorted(lam)]; assert v == sorted(v)
for s in ('problem','answer','solution'):
    assert '400' in it['problem']

print('全検証OK:', ITEM)
