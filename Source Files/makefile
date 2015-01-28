# simple version 1.1
default: clean

gmb:
	gcc -o gmb -O3 gmb.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
	
cfdu:
	gcc -o cfdu -O3 cfdu.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl


listdir:
	gcc -o listdir listdir.c -lm -lnsl


	
gmbdouble:
	gcc -o gmb-double gmb-double.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
compare:
	gcc -o compare compare.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
nmf:
	gcc -o nmf nmf.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
bmf:
	gcc -o bmf bmf.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl

gmf:
	gcc -o gmf gmf.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
rrfx:
	gcc -o rrf rrf.c -I../cfitsio -L../cfitsio -lcfitsio -lm
rrf:
	gcc -o rrf rrf.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	gcc -o rrf rrf.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
rrf2:
	gcc -o rrf2 rrf2.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
rrffast:
	gcc -o rrf -O3 rrf.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl


gmf-test:
	gcc -o gmf-test gmf-test.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl

subtractbias:
	gcc -o subtractbias subtractbias.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
reduction:
	gcc -o reduction reduction.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
pdra:
	gcc -o pdra pdra.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
masterbias:
	gcc -o masterbias masterbias.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl

imarith2:
	gcc -o imarith2 imarith2.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
imarith:
	gcc -o imarith imarith.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
imlist:
	gcc -o imlist imlist.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl


listhead: 
	gcc -o listhead listhead.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
clean:
	rm masterbias
	rm mb.fits
	
