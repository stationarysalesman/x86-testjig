EXECUTABLE = testlib

HDRS = \
	hmmer/src/base/p7_bg.h\
	hmmer/src/base/p7_profile.h\
	hmmer/src/base/p7_hmm.h\
	hmmer/src/base/p7_trace.h\
	hmmer/src/base/p7_anchors.h\
	hmmer/src/base/p7_hmmfile.h\
	hmmer/src/base/p7_anchorhash.h\
	hmmer/src/base/p7_envelopes.h\
	hmmer/src/base/general.h\
	hmmer/src/base/p7_hmmwindow.h\
	hmmer/src/dp_sparse/p7_sparsemx.h\
	hmmer/src/dp_sparse/p7_engine.h\
	hmmer/src/dp_sparse/sparse_viterbi.h\
	hmmer/src/dp_sparse/sparse_fwdback.h\
	hmmer/src/dp_sparse/sparse_decoding.h\
	hmmer/src/dp_sparse/sparse_anchors.h\
	hmmer/src/dp_sparse/sparse_asc_fwdback.h\
	hmmer/src/dp_sparse/sparse_envelopes.h\
	hmmer/src/dp_sparse/sparse_trace.h\
	hmmer/src/dp_sparse/sparse_null2.h\
	hmmer/src/dp_sparse/sparse_aec_align.h\
	hmmer/src/dp_vector/p7_checkptmx.h\
	hmmer/src/dp_vector/msvfilter.h\
	hmmer/src/dp_sparse/p7_spascmx.h\
	hmmer/src/dp_vector/p7_oprofile.h\
	hmmer/src/dp_vector/p7_filtermx.h\
	hmmer/src/dp_vector/simdvec.h\
	hmmer/src/search/modelconfig.h\
	hmmer/src/search/p7_mpas.h\
	hmmer/src/build/modelsample.h\
	hmmer/src/build/modelstats.h\
	hmmer/src/build/evalues.h\
	hmmer/src/misc/emit.h\
	hmmer/src/misc/logsum.h

OBJS = ${HDRS:.h=.o}

CC     = gcc
CFLAGS = -O0 -g
ESLDIR = easel
AR = /usr/bin/ar
MYLIBDIRS = -L./${ESLDIR}
MYSOURCEDIRS= -I./${ESLDIR} -I./hmmer/src/
RANLIB = ranlib
myexe: ${OBJS}
	${CC} ${CFLAGS} ${MYLIBDIRS} ${MYSOURCEDIRS} -o $@ ${OBJS} -lm -leasel  

#all:   px px_serial
#all: ${SOURCES} ${EXECUTABLE}

#${EXECUTABLE}: ${OBJECTS}
#	${CC} ${CFLAGS} ${LDFLAGS} ${OBJS} -o $@

.FORCE:

${OBJS}: ${HDRS} hmmer/src/p7_config.h

#libtest: ${OBJS}
#	${CC} ${CFLAGS} -o libtest -lm -leasel ${OBJS}
.c.o:
	${CC} -c ${CFLAGS} ${MYLIBDIRS} ${MYSOURCEDIRS} -o $@ $< -lm -leasel 

libtest.a: ${OBJS}
	${AR} -r libtest.a ${OBJS} 
	@${RANLIB} libtest.a

#px:     px.c
#	${CC} ${CFLAGS} -o px -L ${HOME}/Documents/research/hmmer-port/code/hmmer/src -L ${HOME}/Documents/research/hmmer-port/code/easel -I ${HOME}/Documents/research/hmmer-port/code/hmmer/src -I ${HOME}/Documents/research/hmmer-port/code/easel px.c -leasel -lm -lpthread

px: px.c
	${CC} ${CFLAGS} -o px px.c ${MYLIBDIRS} -ltest -leasel -lm -lpthread

px_serial: px_serial.c
	${CC} ${CFLAGS} ${MYLIBDIRS} ${MYSOURCEDIRS} -o px_serial px_serial.c -L. -lhmmer -leasel -lm -lpthread

#px_serial:  px_serial.c
#	${CC} ${CFLAGS} -o px_serial -L ${HOME}/Documents/research/hmmer-port/code/hmmer/src -L ${HOME}/Documents/research/hmmer-port/code/easel -I ${HOME}/Documents/research/hmmer-port/code/hmmer/src -I ${HOME}/Documents/research/hmmer-port/code/easel px_serial.c -leasel -lm -lpthread

clean:
	-rm *.o *~
	-rm px px_serial
