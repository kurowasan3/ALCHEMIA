# INO-TRN-045 ヘモグロビンによる酸素の運搬量と溶けた酸素の比較
from sympy import Rational as R
Hb=R(150)/64500; heme=4*Hb
print('Hb mol/L %.4g heme %.4g'%(Hb,heme))
s_l,s_r,s_e=R(97,100),R(70,100),R(30,100)
print('bound lung %.4g'%(heme*s_l))
rest=heme*(s_l-s_r); ex=heme*(s_l-s_e); print('released rest %.4g exercise %.4g ratio %.3g'%(rest,ex,ex/rest))
k=R(1,10**8); d_l=k*13000; d_r=k*5300; d_e=k*2700
print('dissolved lung %.3g, released rest %.3g, ex %.3g'%(d_l,d_l-d_r,d_l-d_e))
print('ratio rest Hb/dissolved %.3g'%(rest/(d_l-d_r)))
print('anemia 100 g/L rest release %.3g ratio %.3g'%(rest*R(100,150),R(100,150)))
