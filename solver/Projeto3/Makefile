CC = g++
DFLAGS = -g -O3
OBJDIR = ../../obj
OBJS = $(OBJDIR)/p3.o $(OBJDIR)/p3src.o \
	$(OBJDIR)/mesh.o $(OBJDIR)/linearsyssolvers.o
IDIR = ../../inc
INCS = $(IDIR)/Mesh.h $(IDIR)/json.hpp $(IDIR)/p3Header.h $(IDIR)/LinearSysSolvers.h
SRC = ../../src
EXEC = p3

linkage: fcheck $(OBJS)
	$(CC) -o $(EXEC) $(OBJS)

fcheck:
	@if ! [ -d "$(OBJDIR)" ]; then mkdir $(OBJDIR); fi

$(OBJDIR)/p3.o: p3.cpp $(INCS)
	$(CC) $(DFLAGS) -I$(IDIR) -c p3.cpp -o $@

$(OBJDIR)/mesh.o: $(SRC)/mesh.cpp $(IDIR)/Mesh.h
	$(CC) $(DFLAGS) -I$(IDIR) -c $(SRC)/mesh.cpp -o $@

$(OBJDIR)/p3src.o: $(SRC)/p3src.cpp $(INCS)
	$(CC) $(DFLAGS) -I$(IDIR) -c $(SRC)/p3src.cpp -o $@

$(OBJDIR)/linearsyssolvers.o: $(SRC)/linearsyssolvers.cpp $(IDIR)/Operators.h $(IDIR)/LinearSysSolvers.h
	$(CC) $(DFLAGS) -I$(IDIR) -c $(SRC)/linearsyssolvers.cpp -o $@

.PHONY: clean
clean:
	-rm $(EXEC)
	-rm $(OBJDIR)/*.o
	-rm -r $(OBJDIR)
