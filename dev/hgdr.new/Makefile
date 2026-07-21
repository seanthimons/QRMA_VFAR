########################################################
# 20060608, AvK:
# * alle namen geprfixed met mth_ (om conflicten met stdlib-functies
#	onmogelijk te maken)
# * idem, preprocessor-constanten: MTH_
# * een extra van preprocessor-#defines WRAP_xxx aangebracht.
#	dit zit in mth_wrap.h
#	je kan nu per functie instellen of je een library-functie,
#	of een user-supplied functie inplugt.
#	Uiteraard geldt dit niet alleen voor de extern zichtbare
#	functies, maar juist ook degenen, die intern gebruikt worden.
#	(fabs, floor, exp). Is voornamelijk bedoeld om te testen.
# * cruft-removal&cleanup
# * driver.c toegevoegd. Doet niet meer dan een minimale
#	functionele test van exp() gamma() en hyperg(), YMMV...
# * Momenteel zijn -strikt genomen- alleen gamma(), lgamma() en hyperg()
#	nodig. De rest zit in voldoende mate in libm.
#	Ik zou toch aanraden om de niet-gebruikt code te laten staan,
#	om bij calamiteiten te kunnen verdeel&heersen.
########################################################
#
SYS=MacOSX-x86
# voor profiling : -pg
# daarna gprof executablefilenaam
# (zit hier in CC-macro, omdat de linker hem ook nodig heeft)
# CC=gcc -std=c99 -pg
CC=gcc -std=c99
# CFLAGS= -Wall -g3 -O6 -pg
CFLAGS= -Wall -g3 -O6
########################################################
# Don't trust the dependencies in this Makefile too much.
# Most of the objects depend on ptm_wrap.h
# (which is included via protos.h, or vice versa ;-) 
# , so after changing ptm_wrap.h, you should
# run 'make clean' , before making.
# BTW: utf8.c is a filter to avoid the silly high ascii
# values, which gcc wants to barf into stderr.
########################################################

MOBJS = hyperg.o gamma.o pow.o powi.o polevl.o \
	floor.o const.o mtherr.o

OBJS = hypg.o ${MOBJS}
#

####
.c.o :
	( ${CC} -c ${CFLAGS} $< ) 2>&1 | ./utf8

###################
all: utf8 hypgdr driver

clean:
	@-rm ${OBJS} hypgdr driver gmon.out

###################
hypgdr: ${OBJS}
	${CC} ${OBJS} -lm -lstdc++ -o $@

driver: driver.o ${MOBJS}
	${CC} driver.o ${MOBJS} -lm -lstdc++ -o $@

################# Eof
