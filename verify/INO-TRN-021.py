# INO-TRN-021 検証：ステンレス鋼（Fe・Cr・Ni）の分析
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
D = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))['results']
# 操作1：過剰 NaOH で Fe(OH)3, Ni(OH)2 は溶けずに残る（Cr は H2O2 で CrO4^2- になってろ液へ）
assert D['Fe^3+']['NaOH-ex']['kind'] == 'remain' and D['Ni^2+']['NaOH-ex']['kind'] == 'remain'
# 操作2：過剰 NH3 で Ni は [Ni(NH3)6]2+（青紫）として溶け，Fe(OH)3（赤褐色）は残る
assert D['Ni^2+']['NH3-ex']['f'] == '[Ni(NH3)6]^2+' and D['Ni^2+']['NH3-ex']['color'] == '青紫'
assert D['Fe^3+']['NH3-ex']['kind'] == 'remain' and D['Fe^3+']['NaOH']['color'] == '赤褐色'
# 操作3：PbCrO4（黄）
assert D['Pb^2+']['CrO4']['f'] == 'PbCrO4' and D['Pb^2+']['CrO4']['color'] == '黄'
# (4) Cu を含む場合：Cu(OH)2 は過剰 NaOH に溶けず P に入り，過剰 NH3 で [Cu(NH3)4]2+ として T へ
assert D['Cu^2+']['NaOH-ex']['kind'] == 'remain' and D['Cu^2+']['NH3-ex']['f'] == '[Cu(NH3)4]^2+'
# 酸性で H2S：CuS は沈殿，Ni2+ は沈殿しない
assert D['Cu^2+']['H2S-a']['f'] == 'CuS' and D['Ni^2+']['H2S-a']['kind'] == 'none'
# 酸化数
from fractions import Fraction
def ox_Cr(ion):  # CrO4^2- : x + 4(-2) = -2
    c, q = parse_formula(ion); return q + 2 * c['O']
assert ox_Cr('CrO4^2-') == 6
for eq in ['2Cr^3+ + 3H2O2 + 10OH- -> 2CrO4^2- + 8H2O', 'Pb^2+ + CrO4^2- -> PbCrO4', 'Ni(OH)2 + 6NH3 -> [Ni(NH3)6]^2+ + 2OH-',
           'Cu(OH)2 + 4NH3 -> [Cu(NH3)4]^2+ + 2OH-', 'Cu^2+ + H2S -> CuS + 2H+', 'Fe^3+ + 3OH- -> Fe(OH)3']:
    assert check_equation(eq) is None, eq
# H2O2 の電子数：Cr 2 mol で 6e-，H2O2 3 mol で 6e-
assert 2 * (6 - 3) == 3 * 2
print('INO-TRN-021 OK')
