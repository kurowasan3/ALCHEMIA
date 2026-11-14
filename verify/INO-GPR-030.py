# 検証スクリプト（自動生成ヘッダ＋問題固有の検証）  python3 -I <this file>
import json, os, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
ID = os.path.splitext(os.path.basename(__file__))[0]
ITEM = json.load(open(f'/home/claude/chem/out/items/{ID}/item.json', encoding='utf-8'))
def _ce(s):
    out = []; i = 0
    while (i := s.find('\\ce{', i)) >= 0:
        j = i + 4; dep = 1
        while dep:
            dep += {'{': 1, '}': -1}.get(s[j], 0); j += 1
        out.append(s[i+4:j-1]); i = j
    return out
EQS = []
for _k in ('problem', 'answer', 'solution'):
    for _c in _ce(ITEM[_k]):
        if len(re.findall(r'<=>|->|<-', _c)) == 1 and not re.search(r'\\boxed|\\phantom|\\text|\bR\b|\bn\b', _c):
            EQS.append(_c)
for _c in EQS:
    _r = check_equation(_c)
    assert _r is None, (_c, _r)
print(f'{ID}: 反応式 {len(EQS)} 本はすべて原子数・電荷がつり合っている')
M = dict(H=1.0, He=4.0, C=12.0, N=14.0, O=16.0, Na=23.0, Mg=24.0, Al=27.0, Si=28.0, P=31.0, S=32.0, Cl=35.5, K=39.0,
         Ca=40.0, Ti=48.0, Cr=52.0, Mn=55.0, Fe=56.0, Co=59.0, Cu=63.5, Zn=65.4, Ga=69.7, Br=80.0, Ag=108.0, I=127.0, Pt=195.0, Au=197.0, Pb=207.0)
def mass(f, masses=None):
    m = masses or M
    c, _ = parse_formula(f)
    return sum(m[e] * float(n) for e, n in c.items())
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
# ---- 問題固有の検証 ----
assert len(EQS) >= 6, EQS
# 加熱の要否・キップの装置（教科書の標準的な整理）を表にして答えと照合
table = {  # 記号: (加熱, 固体のみ, 塊状固体+液体)
    'a': (True, True, False), 'b': (False, False, True), 'c': (True, False, False), 'd': (True, False, False),
    'e': (False, False, True), 'f': (True, False, False), 'g': (True, True, False), 'h': (False, False, True)}
heat = [k for k, v in table.items() if v[0]]
solid = [k for k, v in table.items() if v[1]]
kipp = [k for k, v in table.items() if v[2] and not v[0]]
print(heat, solid, kipp)
ans = ITEM['answer']
assert '(1) ' + '，'.join(heat) in ans and '(2) ' + '，'.join(solid) in ans and '(3) ' + '，'.join(kipp) in ans
# 図鑑の気体の製法で加熱条件を確認
import json
R = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json'))
RX = {r['id']: r for r in R['reactions']}
for i in ['gas-nh3', 'gas-cl2-mno2', 'gas-hcl', 'gas-co', 'gas-o2-kclo3']:
    assert '加熱' in RX[i]['cond'], i
for i in ['gas-h2', 'gas-h2s']:
    assert '加熱' not in RX[i]['cond'], i
co2 = [r for r in R['reactions'] if r.get('gas') == 'CO2' and r.get('cat') == 'gas']
print([ (r['eq'], r['cond']) for r in co2])
assert any('加熱' not in r['cond'] and 'CaCO3' in r['eq'] for r in co2)
print('ALL OK')
