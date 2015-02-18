CC=icc
R_PREFIX=$(shell which R | xargs dirname | xargs dirname)/lib64/R
R_INCLUDE=$(R_PREFIX)/include
R_LIBS=$(R_PREFIX)/lib 

bet : src/bet.o src/geodist.o src/utils.o
	$(CC) -std=gnu99 -shared -L$(R_LIBS) -L/usr/local/lib64 -lR -o lib/bet.so src/bet.o src/geodist.o src/utils.o 	

%.o : %.c
	$(CC) -std=gnu99 -I$(R_INCLUDE) -DNDEBUG  -I/usr/local/include  -fpic  -g -O2  -c $< -o $@

clean :
	rm -f lib/bet.so src/*.o
