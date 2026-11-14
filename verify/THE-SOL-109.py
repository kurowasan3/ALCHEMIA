# THE-SOL-109 verification: octanol/water partition coefficient
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve
P=10**0.88; print('P butanol',P)
frac=P*10/(P*10+100); print('fraction in octanol',frac,'g',1.00*frac)
assert round(frac,2)==0.43
V=symbols('V',positive=True); Vo=solve(P*V/(P*V+100)-0.90,V)[0]; print('Vo for 90%',Vo)
assert round(float(Vo),-1)==120
# two extractions with 59.3 mL? (compare) -- two times 10 mL
rem=(100/(P*10+100))**2; print('two 10 mL extractions remaining',rem,'extracted',1-rem)
BCF5=0.050*10**5; BCF2=0.050*10**2
print('BCF',BCF5,BCF2,'fish ng/kg',2.0*BCF5,2.0*BCF2)
assert BCF5==5000 and 2.0*BCF5==1.0e4 and abs(2.0*BCF2-10)<1e-9
# log P ordering (literature values, for reference): glycerol -1.76, ethanol -0.31, 1-butanol 0.88, toluene 2.73, 1,2,4-trichlorobenzene 4.02
vals={'glycerol':-1.76,'ethanol':-0.31,'1-butanol':0.88,'toluene':2.73,'TCB':4.02}
assert sorted(vals,key=vals.get)==['glycerol','ethanol','1-butanol','toluene','TCB']
assert check_equation('C6H5OH + NaOH -> C6H5ONa + H2O') is None
print('OK')
