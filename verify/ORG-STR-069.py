# ORG-STR-069 正二十面体（カルボラン B10C2H12）の頂点配置異性体
import itertools, numpy as np, networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
p=(1+5**0.5)/2
V=[]
for a in (1,-1):
    for b in (p,-p):
        V+= [(0,a,b),(a,b,0),(b,0,a)]
V=np.array(V,float); n=len(V)
G=nx.Graph()
d=np.linalg.norm(V[:,None]-V[None],axis=2)
for i in range(n):
    for j in range(i+1,n):
        if abs(d[i,j]-2)<1e-6: G.add_edge(i,j)
print('vertices',n,'edges',G.number_of_edges(),'degree',set(dict(G.degree()).values()))
autos=[m for m in GraphMatcher(G,G).isomorphisms_iter()]
def det(m):
    A=np.linalg.lstsq(V,V[[m[i] for i in range(n)]],rcond=None)[0]; return round(np.linalg.det(A))
rot=[m for m in autos if det(m)==1]
print('automorphisms',len(autos),'rotations',len(rot))
sp=dict(nx.shortest_path_length(G,0)); from collections import Counter; print('distance counts from v0',Counter(sp.values()))
def classes(labels_list,group):
    seen=set();cl=[]
    for L in labels_list:
        if L in seen: continue
        orb={tuple(L[m_inv] for m_inv in [inv[i] for i in range(n)]) for inv in [{v:k for k,v in g.items()} for g in group]}
        seen|=orb; cl.append(L)
    return cl
def lab(Cs,Cl=None):
    return tuple('C' if i in Cs else ('X' if i==Cl else 'B') for i in range(n))
two=[lab(set(c)) for c in itertools.combinations(range(n),2)]
cl=classes(two,autos); print('B10C2H12 isomers (full group)',len(cl),'rotation group',len(classes(two,rot)))
for L in cl:
    Cs=[i for i,x in enumerate(L) if x=='C']; dist=nx.shortest_path_length(G,Cs[0],Cs[1])
    mono=[lab(set(Cs),b) for b in range(n) if b not in Cs]
    full=classes(mono,autos); rotc=classes(mono,rot)
    monoC=[lab(set(Cs),Cs[0]),lab(set(Cs),Cs[1])]
    print('C-C distance',dist,'B-Cl isomers',len(full),'(with enantiomers',len(rotc),')','C-Cl isomers',len(classes(monoC,autos)))
# 解説の記述（各組の大きさ・C との隣接数・C の真反対かどうか）の確認
anti={i:int(np.argmin(-d[i])) for i in range(n)}
for c2 in [1,2,3]:
    Cs=(0,[j for j in range(n) if sp[j]==c2][0])
    stab=[g for g in autos if {g[Cs[0]],g[Cs[1]]}==set(Cs)]
    seen=set()
    for b in range(n):
        if b in Cs or b in seen: continue
        orb={g[b] for g in stab}; seen|=orb
        print('dist',c2,'orbit size',len(orb),'adjacent C',sum(G.has_edge(b,c) for c in Cs),'antipode of C',any(anti[c]==b for c in Cs))
