### Intel OneAPI compiler environment ###
F77 = ifort
F77FLAGS += -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3
F90 = ifort
F90FLAGS += -O3 -ip -no-prec-div -qopenmp -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3
#TARGET_CPU = -march=core-avx2
LIB = -L${MKLROOT}/lib/intel64 -qmkl=parallel -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
### ###




### laurel3, cinnamon3 in Kyoto Univ., Intel compiler + IMSL ###
#F77 = ifort
#F90 = ifort
#F90FLAGS += -qopenmp -O3 -ftz -mcmodel=medium -shared-intel
##F90FLAGS += -openmp -check bounds
#TARGET_CPU = -xHost
#LIB = $$LINK_FNL
##LIB = -qmkl
### ###

### ohtaka in ISSP, Intel compiler + Intel MKL ###
### used module ###
#### module load oneapi_compiler/2023.0.0 oneapi_mkl/2023.0.0 oneapi_mpi/2023.0.0
F77 = ifort
F77FLAGS += -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3
F90 = ifort
F90FLAGS += -O3 -ip -no-prec-div -qopenmp -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3
TARGET_CPU = -march=core-avx2
LIB = -L${MKLROOT}/lib/intel64 -qmkl=parallel -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
### ###



PROG := qedynamics_dirac
BIN_DIR := bin
SRC_DIR := src
OBJ_DIR := obj

DIRF = codeF
DIRS = codeS
F90_SRC_FILE =\
	module.f90\
	$(DIRF)/moduleF.f90\
	$(DIRS)/moduleS.f90\
	relvfunc.f90\
	simQED.f90\
	sub_arrangeEig.f90\
	sub_density2.f90\
	sub_density_nuc.f90\
	sub_derivative.f90\
	sub_diffsolver.f90\
	sub_eigsym.f90\
	sub_exvalphys.f90\
	sub_gauss_int.f90\
	sub_int.f90\
	sub_int_nuc.f90\
	sub_intret.f90\
	sub_int_Sch.f90\
	sub_makeQmat.f90\
	sub_norm_MRDFT.f90\
	sub_readBasis.f90\
	sub_readDiracOutput.f90\
	sub_readDiracOutput_KRCI.f90\
	sub_readVectors.f90\
	sub_readVectors_KRCI.f90\
	sub_twoele_divided.f90

F90_DIRF_SRC_FILE =\
	sub_AM.f90\
	sub_Arad.f90\
	sub_calcprop.f90\
	sub_CIreadnatu.f90\
	sub_derivative_AM.f90\
	sub_edmnucchk.f90\
	sub_gauss_int_f.f90\
	sub_int_phys.f90\
	sub_int_Qphys.f90\
	sub_intEDMelefast.f90\
	sub_KRCI3.f90\
	sub_perturbation.f90\
	sub_physqua.f90\
	sub_physqua_2c.f90\
	sub_pt_Qmat.f90\
	sub_save_calEcalC2.f90\
	sub_tdpt.f90

F90_DIRS_SRC_FILE =\
	sub_calc_EDM.f90\
	sub_DIRRCI_EDM.f90\
	sub_integral_EDM.f90\
	sub_local_EDM.f90

F77_SRC_FILE = sub_eig_ch.f

F90_SRC = $(addprefix $(SRC_DIR)/, $(F90_SRC_FILE)) \
		  $(addprefix $(SRC_DIR)/${DIRF}/, $(F90_DIRF_SRC_FILE)) \
		  $(addprefix $(SRC_DIR)/${DIRS}/, $(F90_DIRS_SRC_FILE))
F77_SRC = $(addprefix $(SRC_DIR)/, $(F77_SRC_FILE))

F90_OBJS = $(F90_SRC:src/%.f90=$(OBJ_DIR)/%.o)
F77_OBJS = $(F77_SRC:src/%.f=$(OBJ_DIR)/%.o)
OBJS = $(F90_OBJS) $(F77_OBJS)

MOD_FLAGS = -module $(OBJ_DIR) -I $(OBJ_DIR)

### DO NOT MODIFY THIS BIN_DIR and OBJ_DIR!
### This is a safety measure to prevent variables from becoming empty. ###
BIN_DIR ?= bin
OBJ_DIR ?= obj

.PHONY: all clean show_objs

all: make_dirs $(BIN_DIR)/$(PROG)

make_dirs:
	mkdir -p $(BIN_DIR) $(OBJ_DIR) $(OBJ_DIR)/${DIRF} $(OBJ_DIR)/${DIRS}

show_objs:
	@echo $(OBJS)
	@echo $(F90) $(TARGET_CPU) $(F90FLAGS)

$(BIN_DIR)/$(PROG): $(OBJS)
	$(F90) $(TARGET_CPU) $(F90FLAGS) $^ $(LIB) $(MOD_FLAGS) -o $@

$(OBJ_DIR)/%.o: src/%.f
	$(F77) $(F77FLAGS) $(MOD_FLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: src/%.F
	$(F77) $(F77FLAGS) $(MOD_FLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: src/%.F90
	$(F90) $(F90FLAGS) $(INCS) $(MOD_FLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: src/%.f90
	$(F90) $(F90FLAGS) $(INCS) $(MOD_FLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR)
	rm -rf $(BIN_DIR)
	rm -f *.mod *.log *~
