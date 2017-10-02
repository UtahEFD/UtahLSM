COMP = gfortran

VPATH = src:src/modules:src/interfaces:

OBJ  = 	LSMmodules.o LSM_ShingletonV1.o readInputs.o \
	solveGroundBC.o netSurfaceRadiation.o integrateSoilDiffusion.o tridag.o \
        getSoilThermalTransfer.o getWaterConductivity.o \
        getStabilityCorrections.o getSurfaceMixingRatio.o

# general gnu compiler
LSM:    $(OBJ)
	$(COMP) $(OBJ) -o LSM

trace:	$(OBJ1)
	$(COMP) -trace $(OBJ) -o LSM

.f.o:
	$(COMP) -c $?

clean:
	rm -f *.o *.mod LSM
