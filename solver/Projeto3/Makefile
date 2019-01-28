CC = g++
DFLAGS = -g
OBJDIR = ../../obj
OBJS = $(OBJDIR)/test.o $(OBJDIR)/flowSolver.o $(OBJDIR)/mecflu.o \
	$(OBJDIR)/mesh.o $(OBJDIR)/numerics.o
IDIR = ../../inc
INCS = $(IDIR)/MecFlu.h $(IDIR)/FlowSolver.h $(IDIR)/json.hpp $(IDIR)/Numerics.h
SRC = ../../src
EXEC = test

linkage: fcheck $(OBJS)
	$(CC) -o $(EXEC) $(OBJS)

fcheck:
	@if ! [ -d "$(OBJDIR)" ]; then mkdir $(OBJDIR); fi

$(OBJDIR)/test.o: test.cpp $(INCS)
	$(CC) $(DFLAGS) -I$(IDIR) -c test.cpp -o $@

$(OBJDIR)/flowSolver.o: $(SRC)/flowSolver.cpp $(IDIR)/FlowSolver.h $(IDIR)/MecFlu.h
	$(CC) $(DFLAGS) -I$(IDIR) -c $(SRC)/flowSolver.cpp -o $@

$(OBJDIR)/mecflu.o: $(SRC)/mecflu.cpp $(IDIR)/MecFlu.h
	$(CC) $(DFLAGS) -I$(IDIR) -c $(SRC)/mecflu.cpp -o $@

$(OBJDIR)/mesh.o: $(SRC)/mesh.cpp $(IDIR)/Mesh.h
	$(CC) $(DFLAGS) -I$(IDIR) -c $(SRC)/mesh.cpp -o $@

$(OBJDIR)/numerics.o: $(SRC)/numerics.cpp $(IDIR)/Numerics.h
	$(CC) $(DFLAGS) -I$(IDIR) -c $(SRC)/numerics.cpp -o $@

.PHONY: clean
clean:
	-rm $(EXEC)
	-rm $(OBJDIR)/*.o
	-rm -r $(OBJDIR)
