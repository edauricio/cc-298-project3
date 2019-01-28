CC = g++
DFLAGS = -g
OBJDIR = ../../obj
OBJS = $(OBJDIR)/test.o
IDIR = ../../inc
INCS = $(IDIR)/Bulk.h
SRC = ../../src
EXEC = test

meshgen: fcheck $(OBJS)
	$(CC) -o $(EXEC) $(OBJS)

fcheck:
	@if ! [ -d "$(OBJDIR)" ]; then mkdir $(OBJDIR); fi

$(OBJDIR)/test.o: test.cpp $(INCS)
	$(CC) -I$(IDIR) -c test.cpp -o $@

$(OBJDIR)/aux.o: $(SRC)aux.cpp $(IDIR)Aux.h
	$(CC) -I$(IDIR) -c $(SRC)/aux.cpp -o $@

.PHONY: clean
clean:
	-rm $(EXEC)
	-rm $(OBJDIR)*.o
	-rm -r $(OBJDIR)
