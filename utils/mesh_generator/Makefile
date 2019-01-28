CC = g++
DFLAGS = -g
OBJDIR = ../../obj/
OBJS = $(OBJDIR)mesh_gen.o $(OBJDIR)aux.o $(OBJDIR)mesh.o
IDIR = ../../inc/
INCS = $(IDIR)Mesh.h $(IDIR)Aux.h
SRC = ../../src/
EXEC = meshgen

meshgen: fcheck $(OBJS)
	$(CC) -o $(EXEC) $(OBJS)

fcheck:
	@if ! [ -d "$(OBJDIR)" ]; then mkdir $(OBJDIR); fi

$(OBJDIR)mesh_gen.o: mesh_gen.cpp $(INCS)
	$(CC) -I$(IDIR) -c mesh_gen.cpp -o $@

$(OBJDIR)aux.o: $(SRC)aux.cpp $(IDIR)Aux.h
	$(CC) -I$(IDIR) -c $(SRC)aux.cpp -o $@

$(OBJDIR)mesh.o: $(SRC)mesh.cpp $(IDIR)Mesh.h
	$(CC) -I$(IDIR) -c $(SRC)mesh.cpp -o $@

.PHONY: clean
clean:
	-rm $(EXEC)
	-rm $(OBJDIR)*.o
	-rm -r $(OBJDIR)
