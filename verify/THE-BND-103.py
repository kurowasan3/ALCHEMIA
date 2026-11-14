# THE-BND-103 検証：三方両錐・五方両錐の頂点の2か所を C にした異性体の数，極性，C の配位数
import itertools, math
import numpy as np
def tbp():
    P={1:(0,0,1.2),5:(0,0,-1.2)}
    for k,i in enumerate((2,3,4)):
        a=2*math.pi*k/3; P[i]=(math.cos(a),math.sin(a),0)
    return P
def pbp():
    P={1:(0,0,1.0),7:(0,0,-1.0)}
    for k,i in enumerate((2,3,4,5,6)):
        a=2*math.pi*k/5; P[i]=(math.cos(a),math.sin(a),0)
    return P
def edges(P,tol=0.05):
    d={}
    for a,b in itertools.combinations(P,2):
        d[(a,b)]=np.linalg.norm(np.subtract(P[a],P[b]))
    # 辺＝頂点の隣接：各多面体で決まる組
    return d
def adjacency(P,kind):
    A={i:set() for i in P}
    for a,b in itertools.combinations(P,2):
        za,zb=P[a][2],P[b][2]
        if za!=0 and zb!=0: continue                # 頂点どうし（上下）は結合しない
        if za==0 and zb==0:
            # 赤道の環：隣り合うものだけ（三角形は全部）
            n=len([i for i in P if P[i][2]==0])
            ia=sorted(i for i in P if P[i][2]==0).index(a); ib=sorted(i for i in P if P[i][2]==0).index(b)
            if (ia-ib)%n in (1,n-1): A[a].add(b);A[b].add(a)
        else:
            A[a].add(b);A[b].add(a)                  # 頂点と赤道はすべて結合
    return A
def symmetry_perms(P,A):
    # 隣接関係を保つ置換（グラフの自己同型）＝多面体の対称操作（鏡映を含む）
    V=sorted(P); perms=[]
    for p in itertools.permutations(V):
        m=dict(zip(V,p))
        if all((m[b] in A[m[a]]) for a in V for b in A[a]): perms.append(m)
    return perms
def proper(P,perms):
    # 回転（行列式+1）だけを選ぶ
    V=sorted(P); X=np.array([P[v] for v in V]); out=[]
    for m in perms:
        Y=np.array([P[m[v]] for v in V])
        R,_,_,_=np.linalg.lstsq(X,Y,rcond=None)
        if np.allclose(X@R,Y,atol=1e-6) and np.linalg.det(R)>0: out.append(m)
    return out
for name,P,exp in (('三方両錐',tbp(),3),('五方両錐',pbp(),4)):
    A=adjacency(P,name); perms=symmetry_perms(P,A); rot=proper(P,perms)
    pairs=list(itertools.combinations(sorted(P),2))
    orb_rot=set(); reps=[]
    for pr in pairs:
        key=min(tuple(sorted((m[pr[0]],m[pr[1]]))) for m in rot)
        if key not in orb_rot: orb_rot.add(key); reps.append(pr)
    print(name,'頂点の配位数',{i:len(A[i]) for i in A},'回転の数',len(rot),'異性体',sorted(orb_rot))
    assert len(orb_rot)==exp
    # 鏡像異性体の有無：全対称操作でも同じ数なら鏡像異性体なし
    orb_all=set(min(tuple(sorted((m[a],m[b]))) for m in perms) for a,b in pairs)
    assert len(orb_all)==len(orb_rot)
    for a,b in sorted(orb_rot):
        s=np.add(P[a],P[b]); polar=not np.allclose(s,0)
        print('  C の位置',(a,b),'隣接' if b in A[a] else '非隣接','極性' if polar else '無極性','配位数',len(A[a]),len(A[b]))
print('OK')
