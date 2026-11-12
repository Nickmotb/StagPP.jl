import math
names=["aC","n_fluid","XH2O","XCO2","XCH4","XCO","XH2","nC_fluid","nC_melt","nC_solid","nC_carb","molH_fl","nH_solid",
       "nH2_melt","nH_H2O_melt","a_solid","fH2","fex","S_fluid","it_ideality","vol_fluid","res_ideality"]
rd=lambda f:{int(l.split()[0]):[float(x) for x in l.split()[1:]] for l in open(f)}
J=rd("j_coh.txt"); F=rd("f_coh.txt"); S=rd("states.txt")
def rel(a,b,floor):
    if math.isnan(a) and math.isnan(b): return 0.0
    if math.isnan(a) or math.isnan(b): return float('inf')
    return abs(a-b)/max(abs(a),abs(b),floor)
floors={"aC":1e-30,"fH2":1e-30,"vol_fluid":1e-30,"it_ideality":1,"fex":1e-9,"res_ideality":1}
print(f"{'quantity':12s} {'max rel diff':>12s}  {'state':>6s}   (floor: |x| below ~1e-12 of the budget counts as 0)")
worst={}
for i,n in enumerate(names):
    w=0;ks=0
    for k in J:
        fl=floors.get(n,1e-12)
        r=rel(J[k][i],F[k][i],fl)
        if r>w: w,ks=r,k
    worst[n]=(w,ks); 
    if n!="res_ideality": print(f"{n:12s} {w:12.3e}  {ks:6d}")
itdiff=sum(1 for k in J if J[k][19]!=F[k][19]); print(f"ideality passes differ in {itdiff} of {len(J)} states")
# split by fluid presence for the headline
nfl=sum(1 for k in J if J[k][1]>0); print(f"states with a fluid: {nfl} of {len(J)}")
Jc=rd("j_curve.txt"); Fc=rd("f_curve.txt")
print("C(a) table read, max rel diff at a = 1, 0.1, 0.01, 0.001:",
      [max(rel(Jc[k][i],Fc[k][i],1e-12) for k in Jc) for i in range(4)])
