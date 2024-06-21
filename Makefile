# -
#
# SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
# SPDX-License-Identifier: MIT
#
# -
# Compiler: gnu, intel_cpu, cray_cpu, cray_gpu, nvhpc
ARCH = gnu
# Cases: Boundary Layer, Channel, Taylor-Green vortex, 1D test wave
CASE = -DBL # DBL # DCHA # DTGV
# Equation of state: Ideal Gas (IG), Van der Waals (VdW), Peng-Robinson (PR)
EOS_LAW = -DIG     # DIG # DVdW # DPR
# Transport properties: Constant (IG, VdW, PR), Power Law (IG), Sutherland (IG), JossiStielThodos (VdW), Chung (PR)
VISC_LAW = -DSutherland  # DConstant # DPowerLaw # DSutherland # DJST # DChung
# Benchmark mode: if defined, printing only timesteps without parameters
BENCH =  
# Floating-point numbers
DECOMP_OPTIONS= -DDOUBLE_PREC
# Select FFT type
FFT = fftw3_f03


ifeq ($(ARCH),gnu)
    FLAGS  =  -O0 -march=native # -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -fcheck=all
    FLAGSC  = $(FLAGS) -c -cpp
    COMP = mpif90 -ffixed-line-length-none -std=legacy -J $(MOD)
    LIB =
else ifeq ($(ARCH),intel_cpu)
    FLAGS = -O0 -fpp
    FLAGSC = $(FLAGS) -c
    COMP = mpiifort -132 -module $(MOD)
    LIB =
else ifeq ($(ARCH),cray_cpu)
    FLAGS = -O3 -hfp3 -eZ -hnoacc
    FLAGSC = $(FLAGS) -cc
    COMP = ftn -J  $(MOD)
    LIB =
else ifeq ($(ARCH),cray_gpu)
    FLAGS = -O3 -hfp3 -eZ -hacc -D_AMD
    FLAGSC = $(FLAGS) -cc
    COMP = ftn -D_GPU_DIRECT -J $(MOD)
    LIB =
   # LIB = /users/$$USER/EasyBuild/SW/LUMI-23.03/L/hipfort/0.4-0-cpeCray-23.03-rocm5.2/include/hipfort/amdgcn
else ifeq ($(ARCH),nvhpc)
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

####################################### SIMULATE #######################################

PROGRAM = simulate
SRC = ./src/
OBJ = ./obj/
MOD = ./mod/

OBJS = $(OBJ)factor.o $(OBJ)decomp_2d_constants.o $(OBJ)decomp_2d.o $(OBJ)log.o $(OBJ)glassman.o $(OBJ)io.o  $(OBJ)io_std_units.o \
$(OBJ)nvtx.o $(OBJ)timer.o $(OBJ)math.o \
$(OBJ)param.o $(OBJ)eos_var.o $(OBJ)eos_visc.o $(OBJ)eos.o $(OBJ)comm.o \
$(OBJ)grid.o $(OBJ)finitediff.o $(OBJ)auxl.o $(OBJ)init.o $(OBJ)perturbation.o \
$(OBJ)boundary.o $(OBJ)rhs.o $(OBJ)rk3.o $(OBJ)main_comp.o

all: $(PROGRAM) echo_message

$(PROGRAM): $(OBJS) 
	$(COMP) $(FLAGS) $(OBJS) $(LIB) -o $(PROGRAM)

$(OBJ)main_comp.o: $(SRC)main_comp.f90
	$(COMP)  $(CASE) $(FLAGSC) $(SRC)main_comp.f90 -o $(OBJ)main_comp.o 
$(OBJ)rhs.o: $(SRC)rhs.f90 $(SRC)eulerFlux.f90
	$(COMP)  $(FLAGSC) $(SRC)rhs.f90 -o $(OBJ)rhs.o  
$(OBJ)rk3.o: $(SRC)rk3.f90 
	$(COMP)  $(FLAGSC) $(SRC)rk3.f90 -o $(OBJ)rk3.o
$(OBJ)eos.o: $(SRC)eos.f90 
	$(COMP) $(FLAGSC) $(SRC)eos.f90 -o $(OBJ)eos.o 
$(OBJ)eos_var.o: $(SRC)eos_var.f90
	$(COMP) $(CASE) $(EOS_LAW) $(FLAGSC) $(SRC)eos_var.f90 -o $(OBJ)eos_var.o	
$(OBJ)eos_visc.o: $(SRC)eos_visc.f90
	$(COMP) $(CASE) $(EOS_LAW) $(VISC_LAW) $(FLAGSC) $(SRC)eos_visc.f90 -o $(OBJ)eos_visc.o	
$(OBJ)init.o: $(SRC)init.f90
	$(COMP)  $(EOS_LAW) $(FLAGSC) $(SRC)init.f90 -o $(OBJ)init.o
$(OBJ)grid.o: $(SRC)grid.f90 
	$(COMP)  $(CASE) $(FLAGSC) $(SRC)grid.f90 -o $(OBJ)grid.o
$(OBJ)param.o: $(SRC)param.f90 config_*.h preproc/initBL/inputDNS/initBL_params.h preproc/initCHA/inputDNS/initCHA_params.h
	$(COMP)  $(CASE) $(EOS_LAW) $(VISC_LAW) $(FLAGSC) $(SRC)param.f90 -o $(OBJ)param.o
$(OBJ)finitediff.o: $(SRC)finitediff.f90 
	$(COMP)  $(FLAGSC) $(SRC)finitediff.f90 -o $(OBJ)finitediff.o
$(OBJ)comm.o: $(SRC)comm.f90
	$(COMP)	 $(FLAGSC) $(SRC)comm.f90 -o $(OBJ)comm.o
$(OBJ)boundary.o: $(SRC)boundary.f90
	$(COMP)  $(CASE) $(EOS_LAW) $(FLAGSC) $(SRC)boundary.f90 -o $(OBJ)boundary.o
$(OBJ)perturbation.o: $(SRC)perturbation.f90
	$(COMP)  $(FLAGSC) $(SRC)perturbation.f90 -o $(OBJ)perturbation.o	
$(OBJ)auxl.o: $(SRC)auxl.f90
	$(COMP)  $(CASE) $(FLAGSC) $(SRC)auxl.f90 -o $(OBJ)auxl.o
$(OBJ)math.o: $(SRC)math.f90
	$(COMP)  $(FLAGSC) $(SRC)math.f90 -o $(OBJ)math.o
$(OBJ)io_std_units.o: $(SRC)io_std_units.f90
	$(COMP)  $(OUTPUT_UNIT) $(FLAGSC) $(SRC)io_std_units.f90 -o $(OBJ)io_std_units.o 
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

echo_message:
	@echo "The current architecture is: $(ARCH)"
	@echo "The current EoS is: $(EOS_LAW)"
	@echo "The current TP are: $(VISC_LAW)"
	@echo "The current FP are: $(DECOMP_OPTIONS)"
	@echo "The current FFT is in: $(fftw3_f03)"

####################################### POSTPROCESSING #######################################

OBJS_post = $(OBJ)factor.o $(OBJ)decomp_2d_constants.o $(OBJ)decomp_2d.o $(OBJ)decomp_2d_fft.o $(OBJ)log.o $(OBJ)fft_log.o $(OBJ)glassman.o $(OBJ)io.o $(OBJ)io_std_units.o \
$(OBJ)nvtx.o $(OBJ)timer.o $(OBJ)math.o $(OBJ)param.o $(OBJ)eos_var.o $(OBJ)eos_visc.o $(OBJ)eos.o $(OBJ)comm.o  \
$(OBJ)grid.o $(OBJ)finitediff.o $(OBJ)auxl.o $(OBJ)init.o $(OBJ)perturbation.o \
$(OBJ)boundary.o $(OBJ)rhs.o $(OBJ)rk3.o $(OBJ)auxlpro.o $(OBJ)main_post.o

post: $(OBJS_post)
	$(COMP) $(FLAGS) $(OBJS_post) $(LIB) $(LIBFFT) -o post

$(OBJ)auxlpro.o: $(SRC)post/auxlpro.f90
	$(COMP) $(EOS_LAW) $(FLAGSC) $(SRC)post/auxlpro.f90 -o $(OBJ)auxlpro.o
$(OBJ)fft_log.o: $(SRC)decomp_2d/fft_log.f90
	$(COMP) $(FLAGSC) $(DECOMP_OPTIONS) $(INC) $(SRC)decomp_2d/fft_log.f90 -o $(OBJ)fft_log.o
$(OBJ)decomp_2d_fft.o: $(SRC)decomp_2d/fft_$(FFT).f90
	$(COMP) $(FLAGSC) $(DECOMP_OPTIONS) $(INC) $(SRC)decomp_2d/fft_$(FFT).f90 -o $(OBJ)decomp_2d_fft.o
$(OBJ)main_post.o: $(SRC)post/main_post.f90 $(SRC)param.f90
	$(COMP) $(CASE) $(FLAGSC) $(SRC)post/main_post.f90 -o $(OBJ)main_post.o

####################################### INTERPOLATION #######################################

OBJS_intp = $(OBJ)factor.o $(OBJ)decomp_2d_constants.o $(OBJ)decomp_2d.o $(OBJ)log.o $(OBJ)glassman.o $(OBJ)io.o $(OBJ)io_std_units.o \
$(OBJ)timer.o $(OBJ)math.o \
$(OBJ)param.o $(OBJ)eos_var.o $(OBJ)eos_visc.o $(OBJ)eos.o $(OBJ)comm.o $(OBJ)grid.o $(OBJ)finitediff.o $(OBJ)auxl.o\
$(OBJ)init.o $(OBJ)auxlpol.o $(OBJ)main_interpol.o

interpol: $(OBJS_intp)
	$(COMP) $(FLAGS) $(OBJS_intp) $(LIB) -o interpol
$(OBJ)main_interpol.o: $(SRC)interpol/main_interpol.f90
	$(COMP)  $(CASE) $(FLAGSC) $(SRC)interpol/main_interpol.f90 -o $(OBJ)main_interpol.o
$(OBJ)auxlpol.o: $(SRC)interpol/auxlpol.f90
	$(COMP)  $(CASE) $(FLAGSC) $(SRC)interpol/auxlpol.f90 -o $(OBJ)auxlpol.o

clean:
	$(RM) core mod/* obj/* *.i $(PROGRAM) fort.7 post interpol

cleanfull:
	$(RM) *.i *.txt output/*.txt output/planes/*.* output/stats/*.* output/restart/*.* *.out postproc/results/*.* postproc/results/vort/*.* \
		postproc/results/fft/*.* postproc/results/Yavg/*.* postproc/results/rms/*.* postproc/*.txt \
		postproc/visualize/xplanes/*.* postproc/visualize/yplanes/*.* postproc/visualize/zplanes/*.*

cleanout:
	$(RM) output/*.txt output/planes/*.* output/stats/*.* output/restart/*.* 

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
