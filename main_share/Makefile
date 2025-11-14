CPP     = g++
CCFLAGS = -g -Wall -pedantic -O -I$(HOME)/.local/include
LDFLAGS = -L$(HOME)/.local/lib -lastrocpp -lgsl -lgslcblas -lcfitsio

sources	= \
	dataSample.cpp \
	dataSet.cpp \
	quasarSample.cpp \
	distanceModel.cpp \
	distanceModes.cpp \
	expansionModel.cpp\
	friedmannFit.cpp \
	growthModels.cpp
objects	= $(sources:.cpp=.o)
headers = $(sources:.cpp=.h)

%.o: %.cpp $(headers)
	$(CPP) $(CCFLAGS) -c $< -o $@
%.d: main.x
	./main.x

.PHONY: clean dataclean figclean distclean figs union pantheon

main.x: main.o
	$(CPP) $(CCFLAGS) -o $@ $< $(objects) $(LDFLAGS)
covMatrix.x: covMatrix.o
	$(CPP) $(CCFLAGS) -o $@ $< $(LDFLAGS)
testEdS.x: testEdS.o
	$(CPP) $(CCFLAGS) -o $@ $< $(LDFLAGS)
uncertainty.x: uncertainty.o shearUncertainty.o
	$(CPP) $(CCFLAGS) -o $@ $< shearUncertainty.o $(LDFLAGS)

default: main.x
shear: uncertainty.x
all: default shear

clean:
	$(RM) *.o
auxclean: clean
	$(RM) *.x *.aux *.log *.out data/*.x
distclean: auxclean
	$(RM) *.d *.pdf
	$(RM) -rf fig/union fig/pantheon fig/pantheonBAO fig/*.pdf

dataSample.o: dataSample.h
dataSet.o: dataSet.h dataSample.h
distanceModel.o: distanceModel.h dataSet.o
distanceModes.o: distanceModes.h
expansionModel.o: expansionModel.h
friedmannFit.o: friedmannFit.h
growthModels.o: growthModels.h
main.o: $(objects)
uncertainty.o: shearUncertainty.h
