FC = gfortran
CC = gfortran
FFLAGS = -std=f95

LS = /usr/bin/ls
SRCDIRS  =
SRCEXTS = .f90 .f
SOURCES = $(wildcard $(addprefix *, $(SRCEXTS)))
OBJS    = $(addsuffix .o,$(basename $(SOURCES)))
PROCS = $(basename $(SOURCES))
TARGET = cna2009

.SUFFIXES: .f90 .o .mod

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@
	

$(OBJS): $(SOURCES)	
	$(FC) $(FFLAGS) -c $(SOURCES)
	
	
clean:
	rm -f $(PROCS) $(OBJS)

	