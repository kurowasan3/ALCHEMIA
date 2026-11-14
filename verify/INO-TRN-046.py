# INO-TRN-046 verification: crystal field splitting and colour
hcN=6.63e-34*3.00e8*6.02e23
D=hcN/510e-9/1e3; print('Delta',D); assert round(D)==235
lq=620*(0.194/0.200)**5; print('lambdaQ',lq); assert round(lq)==532
assert 605<=620<=750 and 500<=lq<=560
print('ratio',790/590); assert round(790/590,1)==1.3
d={'Sc3+':21-3-18,'Ti3+':22-3-18,'Ni2+':28-2-18,'Cu2+':29-2-18,'Zn2+':30-2-18,'Ag+':47-1-36}
print(d); assert [k for k,v in d.items() if v in (0,10)]==['Sc3+','Zn2+','Ag+']
print('OK')
