# ORG-HYD-058 ディールス・アルダー反応の生成物と原料（逆合成）・s-シス配座
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
da = AllChem.ReactionFromSmarts('[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1-[C:2]=[C:3]-[C:4]-[C:6]-[C:5]-1')
def DA(diene, dpo):
    out=set()
    for ps in da.RunReactants((Chem.MolFromSmiles(diene), Chem.MolFromSmiles(dpo))):
        p=ps[0]
        try: Chem.SanitizeMol(p); out.add(Chem.MolToSmiles(p))
        except Exception: pass
    return out
targets = {'P 1,2-ジメチルシクロヘキセン':('CC1=C(C)CCCC1','C=C(C)C(C)=C','C=C'),
           'Q 4,5-ジメチルシクロヘキセン':('CC1CC=CCC1C','C=CC=C','CC=CC'),
           'R ノルボルネン':('C1=CC2CCC1C2','C1=CCC=C1','C=C')}
for k,(t,dn,dp) in targets.items():
    got = DA(dn, dp); print(k, CalcMolFormula(Chem.MolFromSmiles(t)), 'from', dn, '+', dp, '->', can(t) in got)
# 誤った組み合わせ（P を butadiene + 2-butene から）では得られないこと
print('butadiene + 2-butene gives P?', can('CC1=C(C)CCCC1') in DA('C=CC=C','CC=CC'))
# シクロペンタジエンの二量化
dimer = DA('C1=CCC=C1','C1=CCC=C1'); print('CPD dimer', {CalcMolFormula(Chem.MolFromSmiles(s)) for s in dimer})
print('CPD', CalcMolFormula(Chem.MolFromSmiles('C1=CCC=C1')))
# 1,4-ペンタジエンは共役していない
m = Chem.MolFromSmiles('C=CCC=C'); print('1,4-pentadiene conjugated?', any(b.GetIsConjugated() for b in m.GetBonds()))
