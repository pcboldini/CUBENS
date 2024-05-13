# -
#
# SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
# SPDX-License-Identifier: MIT
#
# -
# Compiler: gnu, cray_cpu, cray_gpu, nvhpc
TEST = gnu
EOS_LAW = -DIG     # DIG # DVdW # DPR
VISC_LAW = -DSutherland   # DPowerLaw # DSutherland # DJST # DChung
BENCH =  # if defined, printing only timesteps without parameters
DECOMP_OPTIONS= -DDOUBLE_PREC

PROGRAM = simulate
FFT = fftw3_f03

ifeq ($(TEST),gnu)
    FLAGS  =  -O0 -march=native # -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -fcheck=all
    FLAGSC  = $(FLAGS) -c -cpp
    COMP = mpif90 -ffixed-line-length-none -std=legacy -J $(MOD)
    LIB =
else ifeq ($(TEST),cray_cpu)
    FLAGS = -O3 -hfp3 -eP  -hnoacc
    FLAGSC = $(FLAGS) -cc
    COMP = ftn -J  $(MOD)
    LIB =
else ifeq ($(TEST),cray_gpu)
    FLAGS = -O3 -hfp3 -eZ -hacc -D_AMD
    FLAGSC = $(FLAGS) -cc
    COMP = ftn -D_GPU_DIRECT -J $(MOD)
    LIB =
   # LIB = /users/$$USER/EasyBuild/SW/LUMI-23.03/L/hipfort/0.4-0-cpeCray-23.03-rocm5.2/include/hipfort/amdgcn
else ifeq ($(TEST),nvhpc)
    FLAGS = -O0 -fast -acc -cuda -Minfo=accel -D_NVIDIA  #-Mnouniform -Mfprelaxed
    FLAGSC  = $(FLAGS) -c -cpp
    COMP = mpif90 -D_GPU_DIRECT -module $(MOD) #
    LIB = -L/opt/nvidia/hpc_sdk/Linux_x86_64/2023/cuda/lib64/ -lnvToolsExt
endif

ifndef BENCH
        OUTPUT_UNIT=-DOUTPUT_UNIT
else ifdef BENCH
        OUTPUT_UNIT=
endif

# Paths to FFTW 3
ifeq ($(FFT),fftw3_f03)
  FFTW3_PATH = /opt/homebrew/Cellar/fftw/3.3.10_1   # macOS
 # INC=-I$(FFTW3_PATH)/include
  INC=-I/opt/homebrew/Cellar/fftw/3.3.10_1/include
 # LIBFFT=-L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f
  LIBFFT = -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib -lfftw3 -lfftw3f
else ifeq ($(FFT),generic)
  INC=
  LIBFFT=
else ifeq ($(FFT),none)
  INC=
  LIBFFT= 
endif

#######################################
## SIMULATE

SRC = ./src/
OBJ = ./obj/
MOD = ./mod/

OBJS = $(OBJ)factor.o $(OBJ)decomp_2d_constants.o $(OBJ)decomp_2d.o $(OBJ)log.o $(OBJ)glassman.o $(OBJ)io.o  $(OBJ)io_std_units.o \
$(OBJ)nvtx.o $(OBJ)timer.o $(OBJ)math.o \
$(OBJ)param.o $(OBJ)eos_var.o $(OBJ)eos_visc.o $(OBJ)eos.o $(OBJ)table.o $(OBJ)comm.o \
$(OBJ)grid.o $(OBJ)finitediff.o $(OBJ)auxl.o $(OBJ)init.o $(OBJ)perturbation.o \
$(OBJ)boundary.o $(OBJ)rhs.o $(OBJ)rk3.o $(OBJ)main_comp.o

all: $(PROGRAM)

$(PROGRAM): $(OBJS) 
	$(COMP) $(FLAGS) $(OBJS) $(LIB) -o $(PROGRAM)

$(OBJ)main_comp.o: $(SRC)main_comp.f90
	$(COMP)  $(FLAGSC) $(SRC)main_comp.f90 -o $(OBJ)main_comp.o 
$(OBJ)rhs.o: $(SRC)rhs.f90 $(SRC)eulerFlux.f90
	$(COMP)  $(FLAGSC) $(SRC)rhs.f90 -o $(OBJ)rhs.o  
$(OBJ)rk3.o: $(SRC)rk3.f90 
	$(COMP)  $(FLAGSC) $(SRC)rk3.f90 -o $(OBJ)rk3.o
$(OBJ)eos.o: $(SRC)eos.f90 
	$(COMP) $(FLAGSC) $(SRC)eos.f90 -o $(OBJ)eos.o 
$(OBJ)eos_var.o: $(SRC)eos_var.f90
	$(COMP) $(EOS_LAW) $(FLAGSC) $(SRC)eos_var.f90 -o $(OBJ)eos_var.o	
$(OBJ)eos_visc.o: $(SRC)eos_visc.f90
	$(COMP) $(VISC_LAW) $(FLAGSC) $(SRC)eos_visc.f90 -o $(OBJ)eos_visc.o	
$(OBJ)init.o: $(SRC)init.f90
	$(COMP)  $(FLAGSC) $(SRC)init.f90 -o $(OBJ)init.o
$(OBJ)grid.o: $(SRC)grid.f90 
	$(COMP)  $(FLAGSC) $(SRC)grid.f90 -o $(OBJ)grid.o
$(OBJ)param.o: $(SRC)param.f90 config.h initBL/inputDNS/initBL_params.h 
	$(COMP)  $(FLAGSC) -cpp $(SRC)param.f90 -o $(OBJ)param.o
$(OBJ)finitediff.o: $(SRC)finitediff.f90 
	$(COMP)  $(FLAGSC) $(SRC)finitediff.f90 -o $(OBJ)finitediff.o
$(OBJ)comm.o: $(SRC)comm.f90
	$(COMP)	 $(FLAGSC) $(SRC)comm.f90 -o $(OBJ)comm.o
$(OBJ)boundary.o: $(SRC)boundary.f90
	$(COMP)  $(FLAGSC) $(SRC)boundary.f90 -o $(OBJ)boundary.o
$(OBJ)perturbation.o: $(SRC)perturbation.f90
	$(COMP)  $(FLAGSC) $(SRC)perturbation.f90 -o $(OBJ)perturbation.o	
$(OBJ)auxl.o: $(SRC)auxl.f90
	$(COMP)  $(FLAGSC) $(SRC)auxl.f90 -o $(OBJ)auxl.o
$(OBJ)math.o: $(SRC)math.f90
	$(COMP)  $(FLAGSC) $(SRC)math.f90 -o $(OBJ)math.o
$(OBJ)io_std_units.o: $(SRC)io_std_units.f90
	$(COMP)  $(OUTPUT_UNIT) $(FLAGSC) $(SRC)io_std_units.f90 -o $(OBJ)io_std_units.o 
$(OBJ)table.o: $(SRC)table.f90
	$(COMP)  $(FLAGSC) $(SRC)table.f90 -o $(OBJ)table.o	
$(OBJ)nvtx.o: $(SRC)nvtx.f90
	$(COMP)  $(FLAGSC) $(SRC)nvtx.f90 -o $(OBJ)nvtx.o	
$(OBJ)timer.o: $(SRC)timer.f90 $(SRC)nvtx.f90
	$(COMP)  $(FLAGSC) $(SRC)timer.f90 -o $(OBJ)timer.o	
$(OBJ)io.o: $(SRC)decomp_2d/io.f90
	$(COMP) $(FLAGSC) $(DECOMP_OPTIONS) $(SRC)decomp_2d/io.f90 -o $(OBJ)io.o
$(OBJ)factor.o: $(SRC)decomp_2d/factor.f90
	$(COMP) $(FLAGSC) $(DECOMP_OPTIONS) $(SRC)decomp_2d/factor.f90 -o $(OBJ)factor.o
$(OBJ)decomp_2d_constants.o: $(SRC)decomp_2d/decomp_2d_constants.f90
	$(COMP) $(FLAGSC) $(DECOMP_OPTIONS) $(SRC)decomp_2d/decomp_2d_constants.f90 -o $(OBJ)decomp_2d_constants.o
$(OBJ)log.o: $(SRC)decomp_2d/log.f90
	$(COMP) $(FLAGSC) $(DECOMP_OPTIONS) $(SRC)decomp_2d/log.f90 -o $(OBJ)log.o
$(OBJ)decomp_2d.o: $(SRC)decomp_2d/decomp_2d.f90
	$(COMP) $(FLAGSC) $(DECOMP_OPTIONS) $(SRC)decomp_2d/decomp_2d.f90 -o $(OBJ)decomp_2d.o  
$(OBJ)glassman.o: $(SRC)decomp_2d/glassman.f90 
	$(COMP) $(FLAGSC) $(DECOMP_OPTIONS) $(SRC)decomp_2d/glassman.f90 -o $(OBJ)glassman.o


## POSTPROCESSING 

OBJS_post = $(OBJ)factor.o $(OBJ)decomp_2d_constants.o $(OBJ)decomp_2d.o $(OBJ)decomp_2d_fft.o  $(OBJ)log.o $(OBJ)fft_log.o $(OBJ)glassman.o $(OBJ)io.o $(OBJ)io_std_units.o \
$(OBJ)nvtx.o $(OBJ)timer.o $(OBJ)math.o $(OBJ)param.o $(OBJ)eos_var.o $(OBJ)eos_visc.o $(OBJ)table.o $(OBJ)comm.o  \
$(OBJ)grid.o $(OBJ)finitediff.o $(OBJ)auxl.o $(OBJ)eos.o $(OBJ)init.o $(OBJ)perturbation.o \
$(OBJ)boundary.o $(OBJ)rhs.o $(OBJ)rk3.o $(OBJ)postpro.o $(OBJ)main_post.o

post: $(OBJS_post)
	$(COMP) $(FLAGS) $(OBJS_post) $(LIB) $(LIBFFT) -o post

$(OBJ)postpro.o: $(SRC)postpro.f90
	$(COMP) $(EOS_LAW) $(FLAGSC) $(SRC)postpro.f90 -o $(OBJ)postpro.o
$(OBJ)fft_log.o: $(SRC)decomp_2d/fft_log.f90
	$(COMP) $(DECOMP_OPTIONS) -O2 -cpp -c $(INC) $(SRC)decomp_2d/fft_log.f90 -o $(OBJ)fft_log.o
$(OBJ)decomp_2d_fft.o: $(SRC)decomp_2d/fft_$(FFT).f90
	$(COMP) $(DECOMP_OPTIONS) -O2 -cpp -c $(INC) $(SRC)decomp_2d/fft_$(FFT).f90 -o $(OBJ)decomp_2d_fft.o
$(OBJ)main_post.o: $(SRC)main_post.f90 $(SRC)param.f90
	$(COMP)  $(FLAGSC) $(SRC)main_post.f90 -o $(OBJ)main_post.o


## INTERPOLATION TOOL 

OBJS_intp = $(OBJ)factor.o $(OBJ)decomp_2d_constants.o $(OBJ)decomp_2d.o $(OBJ)log.o $(OBJ)io.o $(OBJ)io_std_units.o $(OBJ)timer.o $(OBJ)math.o \
$(OBJ)param.o $(OBJ)eos_var.o $(OBJ)eos_visc.o $(OBJ)comm.o $(OBJ)grid.o $(OBJ)finitediff.o $(OBJ)auxl.o\
$(OBJ)table.o $(OBJ)eos.o $(OBJ)init.o $(OBJ)auxlpol.o $(OBJ)main_interpol.o

interpol: $(OBJS_intp)
	$(COMP) $(FLAGS) $(OBJS_intp) $(LIB) -o interpol
$(OBJ)main_interpol.o: $(SRC)interpol/main_interpol.f90
	$(COMP)  $(FLAGSC) $(SRC)interpol/main_interpol.f90 -o $(OBJ)main_interpol.o
$(OBJ)auxlpol.o: $(SRC)interpol/auxlpol.f90
	$(COMP)  $(FLAGSC) $(SRC)interpol/auxlpol.f90 -o $(OBJ)auxlpol.o


clean:
	$(RM) core mod/* obj/* *.i $(PROGRAM) fort.7 post interpol

cleanfull:
	$(RM) *.i *.txt postproc/planes/*.* restart/*.* *.out postproc/results/*.* postproc/results/vort/*.* postproc/results/fft/*.* postproc/results/Yavg/*.* postproc/*.txt 

cleanpost:
	$(RM) core mod/* obj/* post

cleanpol:
	$(RM) core mod/* obj/* interpol

cleanresu:
	$(RM) postproc/results/*.* postproc/results/fft/*.* postproc/results/vort/*.* postproc/results/Yavg/*.* postproc/results/rms/*.*

cleanvisua:
	$(RM) postproc/visualize/xplanes/*.* postproc/visualize/yplanes/*.* postproc/visualize/zplanes/*.*

cleanfft:
	$(RM) postproc/results/fft/*.*
