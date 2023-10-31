CC = mpicxx
HEADER = bmg.h my_fft.h bmgFlags.h
CFLAGS = -Wall -O3 -std=c++11
INCLUDES = -I/users/PAS1623/yuchiwang/local/fftw/include -I/users/PAS1623/yuchiwang/local/fftw/include
LFLAGS = -L/users/PAS1623/yuchiwang/local/fftw/lib -L/users/PAS1623/yuchiwang/local/fftw/lib
LIBS = -lfftw_mpi -lfftw -lgsl -lgslcblas -lm
SRCS = bmg.c io.c init.c elastic.c evolution.c
OBJS = $(SRCS:.c=.o)
MAIN = bmg
.PHONY: depend clean
#all: $(MAIN)

$(MAIN): $(OBJS) Makefile
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS) 

%.o: %.c $(HEADER) Makefile
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $^
