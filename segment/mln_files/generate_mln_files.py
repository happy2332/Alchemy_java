for i in range(25,801,25):
	f0 = open('segment_mln_int_lifted_'+str(i)+'.txt','w')
	f0.write('#domains\n\
citation={0, ..., '+str(i-1)+'}\n\
position={0, ..., 9}\n\
token={0, ..., '+str(i-1)+'}\n\
field={0, ..., 9}\n\n\
#predicates\n\
Token(token, position, citation)\n\
InField(position, field, citation)\n\
SameField(field, citation, citation)\n\
SameCitation(citation, citation)\n\
Equal(field, field)\n\n\
#formulas\n\
( !Token(t, i, c)   |  InField(i, f, c)) ::1\n\
( !InField(i, f, c) |  InField(i, f, c)) ::3\n\
(  InField(i, f, c) | !InField(i, f, c)) ::2\n\
(  Equal(f1, f2)    | !InField(i, f1, c) | !InField(i, f2, c)) ::5\n\
(  Token(token, position, citation)) ::5\n\
( !InField(position, field, citation)) ::4\n\
( !Equal(field, field)) ::1\n')
	f0.close();
