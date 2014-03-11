include config.mak

BIN=vmc_1 test_1
SRC=Amplitude_1.cpp Timer.cpp RanGen.cpp FileManager.cpp WaveFunction_1.cpp linalg.cpp MetroMC_1.cpp Quantity_1.cpp StagFluxWaveFunction_1.cpp StagFluxTransExciton_1.cpp ArgParse.cpp LatticeStepper_1.cpp StagFluxLongExciton_1.cpp State_1.cpp LatticeState_1.cpp SquareLattice.cpp VectorQuantity_1.cpp ScalarQuantity_1.cpp MatrixQuantity_1.cpp StatSpinStruct_1.cpp ProjHeis_1.cpp defs.cpp
HDR=MetroMC_1.h Amplitude_1.h Quantity_1.h linalg.h RanGen.h Timer.h WaveFunction_1.h FileManager.h BigComplex.h BigDouble.h StagFluxWaveFunction_1.h StagFluxGroundState_1.h StagFluxTransExciton_1.h ArgParse.h LatticeStepper_1.h StagFluxLongExciton_1.h State_1.h LatticeState_1.h Lattice.h SquareLattice.h OverlapTrack_1.h StagMagnTrack_1.h VectorQuantity_1.h ScalarQuantity_1.h MatrixQuantity_1.h StatSpinStruct_1.h ProjHeis_1.h
OBJ=$(SRC:.cpp=.o) zdotu_sub.o

all: $(BIN)

vmc.o: $(HDR) gitversion.h
vmc_1.o: $(HDR) gitversion.h
test_1.o: $(HDR)
test.o: $(HDR)
zdotu_sub.o: zdotu_sub.f
	$(OFORT) -c -o $@ $< -O3
defs.o: defs.h
ProjHeis_1.o: ProjHeis_1.h Quantity_1.h VectorQuantity_1.h MatrixQuantity_1.h LatticeState_1.h Timer.h WaveFunction_1.h Stepper_1.h Amplitude_1.h BigComplex.h BigDouble.h
ScalarQuantity_1.o: Quantity_1.h Timer.h
VectorQuantity_1.o: Quantity_1.h Timer.h
MatrixQuantity_1.o: Quantity_1.h VectorQuantity_1.h Timer.h
StatSpinStruct_1.o: Quantity_1.h VectorQuantity_1.h MatrixQuantity_1.h Stepper_1.h Amplitude_1.h LatticeState_1.h Lattice.h
linalg.o: linalg.h BigComplex.h BigDouble.h blas_lapack.h
RanGen.o: RanGen.h Timer.h
StagFluxWaveFunction_1.o: StagFluxWaveFunction_1.h WaveFunction_1.h linalg.h
StagFluxTransExciton_1.o: StagFluxTransExciton_1.h StagFluxWaveFunction_1.h linalg.h
StagFluxLongExciton_1.o: StagFluxLongExciton_1.h StagFluxWaveFunction_1.h linalg.h
ArgParse.o: ArgParse.h
SquareLattice.o: Lattice.h

.PHONY: gitversion.h
gitversion.h: gitversion.h.tpl
	sh check_gitversion.sh

vmc_1: vmc_1.o $(OBJ)
	$(OCXX) -o $@ $^ $(LDFLAGS)

test_1: test_1.o $(OBJ)
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
