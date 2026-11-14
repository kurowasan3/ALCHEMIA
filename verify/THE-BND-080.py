# --- 反応式の機械的チェック（item.json 中の \ce{...} で矢印が1つの完全な式） ---
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def ce_blocks(s):
    out = []
    i = 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0:
            return out
        k = j + 4; depth = 1
        while depth:
            if s[k] == '{': depth += 1
            elif s[k] == '}': depth -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k

def check_item(path):
    d = json.load(open(path, encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(d[key]):
            if '\\text' in b or '\\mathrm' in b or '\\boxed' in b:
                continue
            arrows = re.findall(r'<=>|->|<-', b)
            if len(arrows) != 1:
                continue
            e = re.sub(r'\[[^\]]*\]', '', b)          # 矢印ラベルを除く
            e = e.replace('<=>', '->')
            e = re.sub(r'\((g|l|s|aq)\)', '', e)
            r = check_equation(e)
            assert r is None, (b, r)
            n += 1
    print('反応式チェック OK:', n, '本')

import itertools
pi, s2, s3 = 3.14, 1.41, 1.73
rp, rm = 0.060, 0.184
# (1) 数え上げ（分率座標で）
corner = 8*(1/8); face = 6*(1/2); assert corner + face == 4
a = 4*(rp+rm)/s3; a3 = round(a, 3); print('a', a, '->', a3); assert a3 == 0.564
dSS = a3/s2; print('S-S', dSS, 'diam', 2*rm); assert round(dSS, 3) == 0.400 and dSS > 2*rm
V = a3**3
f = 4*(4/3)*pi*(rp**3+rm**3)/V; fS = 4*(4/3)*pi*rm**3/V
print('f', f, 'fS', fS); assert round(f*100, 1) == 60.2 and round(fS*100, 1) == 58.1
fd = s3*pi/16; print('diamond', fd); assert round(fd*100, 1) == 34.0
# 幾何の確認（厳密値）：Zn(1/4,1/4,1/4) と S(0,0,0) の距離 = √3a/4，四面体の4頂点まで等距離
import math
zn = (0.25, 0.25, 0.25)
for S in [(0,0,0),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5)]:
    assert abs(math.dist(zn, S) - math.sqrt(3)/4) < 1e-12
check_item('/home/claude/chemx/tokodai/out/items/THE-BND-080/item.json')
print('OK')
