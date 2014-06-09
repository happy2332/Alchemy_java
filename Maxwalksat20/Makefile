maxwalksat: maxwalksat.c
	gcc -O2 maxwalksat.c -lm -o maxwalksat

dist:
	tar cevof maxwalksat.dist  maxwalksat.c README.maxwalksat Makefile
	compress maxwalksat.dist
	uuencode maxwalksat.dist.Z maxwalksat.dist.Z > maxwalksat.dist.Z.uu
	
try:	try.c
	cc try.c -lm -o try

costwalksat: costwalksat.c
	cc costwalksat.c -lm -o costwalksat
