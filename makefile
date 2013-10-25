include config.mak

BIN=test transsqw groundstate longsqw vmc
SRC=SpinState.cpp Amplitude.cpp SpinOp.cpp SpinDensity.cpp Timer.cpp RanGen.cpp FileManager.cpp WaveFunction.cpp linalg.cpp MetroMC.cpp Quantity.cpp StagFluxWaveFunction.cpp StagFluxTransExciton.cpp Jastrow.cpp ArgParse.cpp SpinSpinCorr.cpp StagJastrow.cpp StaggMagnJastrow.cpp ScalarQuantity.cpp VectorQuantity.cpp MatrixQuantity.cpp FullSpaceStepper.cpp ProjHeis.cpp StagFluxLongExciton.cpp StatSpinStruct.cpp
HDR=MetroMC.h SpinState.h Amplitude.h Quantity.h ScalarQuantity.h linalg.h RanGen.h Timer.h WaveFunction.h VectorQuantity.h SpinOp.h MatrixQuantity.h SpinDensity.h FileManager.h Stepper.h BigComplex.h BigDouble.h StagFluxWaveFunction.h StagFluxGroundState.h StagFluxTransExciton.h Jastrow.h ArgParse.h SpinSpinCorr.h StagJastrow.h StaggMagnJastrow.h FullSpaceStepper.h ProjHeis.h StagFluxLongExciton.h StatSpinStruct.h StagMagn.h
OBJ=$(SRC:.cpp=.o) zdotu_sub.o

all: $(BIN)

groundstate.o: $(HDR)
transsqw.o: $(HDR)
test.o: $(HDR)
longsqw.o: $(HDR)
zdotu_sub.o: zdotu_sub.f
	$(OFORT) -c -o $@ $< -O3
MetroMC.o: RanGen.h SpinState.h Stepper.h Quantity.h BigDouble.h Timer.h
SpinState.o: SpinState.h Amplitude.h RanGen.h Timer.h BigComplex.h BigDouble.h
Amplitude.o: Amplitude.h SpinState.h linalg.h Timer.h WaveFunction.h BigComplex.h BigDouble.h blas_lapack.h
WaveFunction.o: WaveFunction.h linalg.h BigComplex.h BigDouble.h Amplitude.h
FullSpaceStepper.o: FullSpaceStepper.h Amplitude.h WaveFunction.h linalg.h RanGen.h Timer.h SpinState.h Stepper.h BigComplex.h BigDouble.h
ProjHeis.o: ProjHeis.h Quantity.h VectorQuantity.h MatrixQuantity.h SpinState.h Timer.h WaveFunction.h Stepper.h Amplitude.h BigComplex.h BigDouble.h
SpinOp.o: SpinState.h Timer.h BigComplex.h BigDouble.h
ScalarQuantity.o: Quantity.h Timer.h
VectorQuantity.o: Quantity.h Timer.h
MatrixQuantity.o: Quantity.h VectorQuantity.h Timer.h
SpinDensity.o: MatrixQuantity.h VectorQuantity.h Quantity.h SpinState.h SpinOp.h Timer.h BigComplex.h BigDouble.h
linalg.o: linalg.h BigComplex.h BigDouble.h blas_lapack.h
RanGen.o: RanGen.h Timer.h
StagFluxWaveFunction.o: StagFluxWaveFunction.h WaveFunction.h linalg.h
StagFluxTransExciton.o: StagFluxTransExciton.h StagFluxWaveFunction.h linalg.h
StagFluxLongExciton.o: StagFluxLongExciton.h StagFluxWaveFunction.h linalg.h
Jastrow.o: Jastrow.h SpinState.h
StagJastrow.o: StagJastrow.h Jastrow.h linalg.h
ArgParse.o: ArgParse.h
SpinSpinCorr.o: SpinSpinCorr.h VectorQuantity.h Quantity.h Stepper.h
StaggMagnJastrow.o: StaggMagnJastrow.h Jastrow.h linalg.h SpinState.h

groundstate: groundstate.o $(OBJ)
	$(OCXX) -o $@ $^ $(LDFLAGS)

test: test.o $(OBJ)
	$(OCXX) -o $@ $^ $(LDFLAGS)

transsqw: transsqw.o $(OBJ)
	$(OCXX) -o $@ $^ $(LDFLAGS)

vmc: vmc.o $(OBJ)
	$(OCXX) -o $@ $^ $(LDFLAGS)

longsqw: longsqw.o $(OBJ)
	$(OCXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp makefile config.mak local.mak
	$(OCXX) -c -o $@ $< $(CFLAGS)

clean:
	rm -f *.o

cleaner: clean
	rm -f $(BIN)

install: $(BIN)
	mkdir -p $(INSTALL_DIR); cp $(BIN) $(INSTALL_DIR); cd collect; make install

showme:
	mpicxx -showme

doc: $(HDR) Doxyfile
	doxygen
