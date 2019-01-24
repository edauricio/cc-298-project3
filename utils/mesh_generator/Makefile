CC = g++
DFLAGS = -g
OBJDIR = ../../obj
OBJS = $(OBJDIR)/mesh_gen.o
IDIR = ../../inc
INCS = $(IDIR)/Mesh.h
EXEC = meshgen

meshgen: $(OBJS)
	$(CC) -o $(EXEC) $(OBJS)

$(OBJDIR)/mesh_gen.o: mesh_gen.cpp $(INCS)
	$(CC) -I$(IDIR) -c mesh_gen.cpp -o $@

.PHONY: clean
clean:
	-rm $(EXEC) $(OBJS)
