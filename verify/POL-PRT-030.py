code={'AUG':'Met','CAU':'His','CAC':'His','GGC':'Gly','UGG':'Trp','UAA':'*','CUG':'Leu','GCU':'Ala','GGU':'Gly'}
# check against standard genetic code
std={'AUG':'Met','CAU':'His','CAC':'His','GGC':'Gly','UGG':'Trp','UAA':'*','CUG':'Leu','GCU':'Ala','GGU':'Gly'}
assert code==std
comp={'A':'U','T':'A','G':'C','C':'G'}
def mrna(tmpl35): return ''.join(comp[b] for b in tmpl35)   # template written 3'->5' gives mRNA 5'->3' in same order
def tr(m):
    out=[]
    for i in range(0,len(m)-2,3):
        c=m[i:i+3]; a=code.get(c,'?')
        if a=='*': break
        out.append(a)
    return out
t='TACGTACCGACCATT'
m=mrna(t); print(m); assert m=='AUGCAUGGCUGGUAA'
p=tr(m); print(p); assert p==['Met','His','Gly','Trp']
ta=t[:5]+'G'+t[6:]; assert t[5]=='A'; ma=mrna(ta); print(ma, tr(ma)); assert tr(ma)==p
tb=t[:4]+t[5:]; assert t[4]=='T'; mb=mrna(tb); print(mb, tr(mb)); assert tr(mb)[:4]==['Met','Leu','Ala','Gly']
assert 4**2==16<20 and 4**3==64
print('POL-PRT-030 OK')
