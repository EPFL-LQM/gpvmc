include config.mak

BIN=vmc test
SRC=Amplitude.cpp Timer.cpp RanGen.cpp FileManager.cpp WaveFunction.cpp linalg.cpp MetroMC.cpp Quantity.cpp StagFluxWaveFunction.cpp StagFluxTransExciton.cpp ArgParse.cpp LatticeStepper.cpp StagFluxLongExciton.cpp State.cpp LatticeState.cpp SquareLattice.cpp VectorQuantity.cpp ScalarQuantity.cpp MatrixQuantity.cpp StatSpinStruct.cpp ProjHeis.cpp defs.cpp StagMagn.cpp StagMagnTrack.cpp
HDR=MetroMC.h Amplitude.h Quantity.h linalg.h RanGen.h Timer.h WaveFunction.h FileManager.h BigComplex.h BigDouble.h StagFluxWaveFunction.h StagFluxGroundState.h StagFluxTransExciton.h ArgParse.h LatticeStepper.h StagFluxLongExciton.h State.h LatticeState.h Lattice.h SquareLattice.h OverlapTrack.h StagMagnTrack.h VectorQuantity.h ScalarQuantity.h MatrixQuantity.h StatSpinStruct.h ProjHeis.h
OBJ=$(SRC:.cpp=.o) zdotu_sub.o

all: $(BIN)

vmc.o: $(HDR) gitversion.h
vmc.o: $(HDR) gitversion.h
test.o: $(HDR)
test.o: $(HDR)
zdotu_sub.o: zdotu_sub.f
	$(OFORT) -c -o $@ $< -O3
defs.o: defs.h
ProjHeis.o: ProjHeis.h Quantity.h VectorQuantity.h MatrixQuantity.h LatticeState.h Timer.h WaveFunction.h Stepper.h Amplitude.h BigComplex.h BigDouble.h
ScalarQuantity.o: Quantity.h Timer.h
VectorQuantity.o: Quantity.h Timer.h
MatrixQuantity.o: Quantity.h VectorQuantity.h Timer.h
StatSpinStruct.o: Quantity.h VectorQuantity.h MatrixQuantity.h Stepper.h Amplitude.h LatticeState.h Lattice.h
linalg.o: linalg.h BigComplex.h BigDouble.h blas_lapack.h
RanGen.o: RanGen.h Timer.h
StagFluxWaveFunction.o: StagFluxWaveFunction.h WaveFunction.h linalg.h
StagFluxTransExciton.o: StagFluxTransExciton.h StagFluxWaveFunction.h linalg.h
StagFluxLongExciton.o: StagFluxLongExciton.h StagFluxWaveFunction.h linalg.h
ArgParse.o: ArgParse.h
SquareLattice.o: Lattice.h

.PHONY: gitversion.h
gitversion.h: gitversion.h.tpl
	sh check_gitversion.sh

vmc: vmc.o $(OBJ)
	$(OCXX) -o $@ $^ $(LDFLAGS)

test: test.o $(OBJ)
	$(OCXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp makefile config.mak local.mak
	$(OCXX) -c -o $@ $< $(CFLAGS)

clean:
	rm -f *.o

cleaner: clean
	rm -f $(BIN)

install: $(BIN) wrap_vmc_debug.sh args_example.in
	mkdir -p $(INSTALL_DIR); cp $(BIN) wrap_vmc_debug.sh args_example.in $(INSTALL_DIR); cd collect; make install

showme:
	mpicxx -showme

doc: $(HDR) Doxyfile
	doxygen
