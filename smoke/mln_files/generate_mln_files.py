for i in range(825,1001,25):
	f0 = open('smoke_mln_'+str(i)+'.txt','w')
	f0.write('#domains\n\
dom1={0, ..., '+str(i-1)+'}\n\n\
#predicates\n\
C(dom1)\n\
S(dom1)\n\
F(dom1,dom1)\n\n\
#formulas\n\
(S(x)) :: 2\n\
(!C(x)) :: 4\n\
(F(x,y)) :: 2\n\
( !S(x) | C(x)) :: 3\n\
( !S(x) | !F(x,y) | S(y)) :: 4')
	f0.close();
