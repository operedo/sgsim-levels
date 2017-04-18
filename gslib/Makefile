FC=$(COMPILER)
FFLAGS=$(FLAGS)
OPENMP=$(OMP)

#FC=gfortran
#FFLAGS= -Ofast  -funroll-loops --param max-unroll-times=4 -ftree-vectorize 
#OPENMP=-fopenmp

#FC=ifort
#FFLAGS=-pg -O3 -march=native
LIB=gslib.a

.SUFFIXES: .o .for

%.o: %.for
	$(FC) $(FFLAGS) $(OPENMP) -c $<

all: $(LIB)

OBJS= \
	acorni.o \
	backtr.o \
	beyond.o \
	blue.o \
	chknam.o \
	chktitle.o \
	cova3.o \
	dlocate.o \
	dpowint.o \
	dsortem.o \
	gauinv.o \
	gcum.o \
	getindx.o \
	getz.o \
	green.o \
	hexa.o \
	ksol.o \
	ktsol.o \
	locate.o \
	nscore.o \
	numtext.o \
	ordrel.o \
	picksupr.o \
	powint.o \
	psfill.o \
	psline.o \
	pstext.o \
	rand.o \
	red.o \
	resc.o \
	scal.o \
	setrot.o \
	setsupr.o \
	sortem.o \
	sqdist.o \
	srchsupr.o \
	strlen.o \

$(LIB): $(OBJS)
	ar rcs $(LIB) $(OBJS)

clean:
	rm $(LIB) $(OBJS)
