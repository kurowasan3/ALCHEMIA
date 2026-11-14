# ORG-ARO-092 C–H のヒドロキシ化とヘミアミナールの分解 検証
from rdkit import Chem
from collections import Counter
def c(s): return Counter(a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(s)).GetAtoms())
def M(s):
    A={'C':12,'H':1.0,'N':14,'O':16}; return sum(A[k]*v for k,v in c(s).items())
S='CN(CC)Cc1ccccc1'; print('S',M(S)); assert M(S)==149
O=Counter({'O':1})
cases={'a N-CH3':('OCN(CC)Cc1ccccc1',['C=O','CCNCc1ccccc1']),
       'b N-CH2':('CN(C(O)C)Cc1ccccc1',['CC=O','CNCc1ccccc1']),
       'c benzyl':('CN(CC)C(O)c1ccccc1',['O=Cc1ccccc1','CNCC']),
       'd CH2OH':('CN(CCO)Cc1ccccc1',[]),
       'e para-OH':('CN(CC)Cc1ccc(O)cc1',[])}
for k,(h,prods) in cases.items():
    assert c(h)==c(S)+O, k
    print(k,'hydroxylated MW',M(h), [ (p,M(p)) for p in prods])
    if prods: assert c(h)==sum((c(p) for p in prods),Counter()), k
assert M('CNCc1ccccc1')==121 and M('CCNCc1ccccc1')==135 and M('CNCC')==59
# アニソールの O-脱メチル化
assert c('COc1ccccc1')+O==c('OCOc1ccccc1') and c('OCOc1ccccc1')==c('C=O')+c('Oc1ccccc1')
print('OK')
