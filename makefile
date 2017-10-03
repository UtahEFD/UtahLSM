# --------------------------------------
# Makefile for UtahLSM
#
# Author: Jeremy Gibbs
# Date: 2017-10-02
#
# UtahLSM developed by Shingleton 2012
# --------------------------------------

TARGET  = LSM
FC      = gfortran
SRCDIR  = src

SOURCES := $(wildcard $(SRCDIR)/*f)
OBJECTS := $(SOURCES:$(SRCDIR)/%.f=$(SRCDIR)/%.o)
OBJECTS := $(SRCDIR)/LSMmodules.o $(SRCDIR)/LSM_ShingletonV1.o \
           $(SRCDIR)/readInputs.o $(SRCDIR)/solveGroundBC.o \
           $(SRCDIR)/netSurfaceRadiation.o $(SRCDIR)/integrateSoilDiffusion.o \
           $(SRCDIR)/tridag.o $(SRCDIR)/getSoilThermalTransfer.o \
           $(SRCDIR)/getWaterConductivity.o $(SRCDIR)/getStabilityCorrections.o \
           $(SRCDIR)/getSurfaceMixingRatio.o

$(TARGET): $(OBJECTS)
	@$(FC) $(OBJECTS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(SRCDIR)/%.o : $(SRCDIR)/%.f
	@$(FC) -J$(SRCDIR) -c $< -o $@
	@echo "Compiled "$<" successfully!"

clean:
	rm -rf $(OBJECTS) $(SRCDIR)/*mod $(TARGET)
	@echo "Cleanup complete!" 