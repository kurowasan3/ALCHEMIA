# ORG-HYD-042 検証：二環式炭化水素の水素化分解（どの C–C 結合が切れるか）で生じる生成物
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as D
RDLogger.DisableLog('rdApp.*')
def cleave(smi):
    m = Chem.MolFromSmiles(smi); res = {}
    for b in m.GetBonds():
        if b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C' and b.IsInRing():
            rw = Chem.RWMol(m); rw.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            p = rw.GetMol(); Chem.SanitizeMol(p); s = Chem.MolToSmiles(p)
            ring = sorted(len(r) for r in p.GetRingInfo().AtomRings())
            res.setdefault(s, []).append((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), ring))
    return res
names = {'C1CCCCC1': 'シクロヘキサン', 'CC1CCCC1': 'メチルシクロペンタン', 'CCC1CC1C': '1-エチル-2-メチルシクロプロパン',
         'C1CCC1': 'シクロブタン', 'CC1CC1': 'メチルシクロプロパン', 'C1CCCC1': 'シクロペンタン', 'CC1CCC1': 'メチルシクロブタン', 'CCCC1CC1': 'プロピルシクロプロパン'}
for smi, nm in [('C1CC2CC2C1', 'ビシクロ[3.1.0]ヘキサン'), ('C1C2CC12', 'ビシクロ[1.1.0]ブタン'), ('C1CC2CC12', 'ビシクロ[2.1.0]ペンタン')]:
    m = Chem.MolFromSmiles(smi); print(nm, D.CalcMolFormula(m))
    r = cleave(smi)
    for s, v in r.items(): print('   ', s, names.get(s, '?'), D.CalcMolFormula(Chem.AddHs(Chem.MolFromSmiles(s))), '切れる結合の数', len(v), '残る環', v[0][2])
    if nm.startswith('ビシクロ[3.1.0]'):
        assert set(r) == {'C1CCCCC1', 'CC1CCCC1', 'CCC1CC1C', 'CCCC1CC1'} and len(r['C1CCCCC1']) == 1 and len(r['CC1CCCC1']) == 2 and len(r['CCC1CC1C']) == 2
    if nm.startswith('ビシクロ[1.1.0]'):
        assert set(r) == {'C1CCC1', 'CC1CC1'}
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['C6H10 + H2 -> C6H12', 'C4H6 + H2 -> C4H8']: assert check_equation(e) is None
print('OK')
