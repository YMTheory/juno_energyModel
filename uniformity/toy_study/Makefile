cc=g++
LIBS=$(shell root-config --libs)
INCLUDE=$(shell root-config --cflags)
exe=nPEMap_corr1
obj=nPEMap_corr1.o
$(exe):$(obj)
	$(cc) $(INCLUDE) -o $(exe) nPEMap_corr1.cc $(LIBS)
$(obj):nPEMap_corr1.cc
	$(cc) -c $(INCLUDE) nPEMap_corr1.cc
clean:
	rm -rf *.o $(exe)
