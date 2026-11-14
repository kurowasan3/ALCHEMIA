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

# --- 数値の検証 ---
from sympy import Rational as R
cCH4, cCO, cH2, vap = -890, -283, -286, 44
dHA = cCH4 - (2*cCO + 2*cH2)
dHB_g = cCO - cH2 - vap
dHB_l = cCO - cH2
print('ΔH_A =', dHA, ' ΔH_B(g) =', dHB_g, ' ΔH_B(l) =', dHB_l)
assert dHA == 248 and dHB_g == -41 and dHB_l == 3
# 酸化数の収支：CH4 の C -4→+2 (+6)，CO2 の C +4→+2 (-2)，H 4個 +1→0 (-4)
assert 6 == 2 + 4
# 気体分子数
assert (1+1, 2+2) == (2, 4) and (1+1, 1+1) == (2, 2)
check_item('/home/claude/chemx/tokodai/out/items/THE-EQL-123/item.json')
print('OK')
