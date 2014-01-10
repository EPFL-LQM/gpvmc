include config.mak

BIN=vmc
SRC=SpinState.cpp Amplitude.cpp SpinOp.cpp SpinDensity.cpp Timer.cpp RanGen.cpp FileManager.cpp WaveFunction.cpp WaveFunction_1.cpp linalg.cpp MetroMC.cpp Quantity.cpp StagFluxWaveFunction.cpp StagFluxWaveFunction_1.cpp StagFluxTransExciton.cpp StagFluxTransExciton_1.cpp Jastrow.cpp ArgParse.cpp SpinSpinCorr.cpp StagJastrow.cpp StaggMagnJastrow.cpp ScalarQuantity.cpp VectorQuantity.cpp MatrixQuantity.cpp FullSpaceStepper.cpp ProjHeis.cpp StagFluxLongExciton.cpp StagFluxLongExciton_1.cpp StatSpinStruct.cpp State_1.cpp
HDR=MetroMC.h SpinState.h Amplitude.h Quantity.h ScalarQuantity.h linalg.h RanGen.h Timer.h WaveFunction.h WaveFunction_1.h VectorQuantity.h SpinOp.h MatrixQuantity.h SpinDensity.h FileManager.h Stepper.h BigComplex.h BigDouble.h StagFluxWaveFunction.h StagFluxWaveFunction_1.h StagFluxGroundState.h StagFluxTransExciton.h StagFluxTransExciton_1.h Jastrow.h ArgParse.h SpinSpinCorr.h StagJastrow.h StaggMagnJastrow.h FullSpaceStepper.h ProjHeis.h StagFluxLongExciton.h StagFluxLongExciton_1.h StatSpinStruct.h StagMagn.h State_1.h
OBJ=$(SRC:.cpp=.o) zdotu_sub.o

all: $(BIN)

vmc.o: $(HDR) gitversion.h
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

.PHONY: gitversion.h
gitversion.h:
	sed "s/\".*\"/\"`git rev-parse --short HEAD`\"/" <gitversion.h >gitversion_temp.h && mv gitversion_temp.h gitversion.h

vmc: vmc.o $(OBJ)
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
