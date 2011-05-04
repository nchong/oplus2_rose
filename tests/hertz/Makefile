COMMONDIR := ../common
OP2ROSE := ../../translator/op2rose

NVCC := nvcc
INCLUDEFLAGS := -I$(COMMONDIR) -I$(HERTZ_COMMONDIR)
NVCCFLAGS := -O2 -arch sm_13 $(INCLUDEFLAGS) -D OP_DIAGS=0 #-g

VPATH = $(COMMONDIR)

SOURCES = rose_hertz.cpp op_par.cu cutil.cpp \
	cmd_arg_reader.cpp stopwatch.cpp bank_checker.cpp stopwatch_linux.cpp $(HERTZ_COMMONDIR)/libcommon.so \

all: src2src build

build: $(SOURCES)
	$(NVCC) $(NVCCFLAGS) $^ -o op2_cpu_timer
	cp rose_hertz.cpp rose_hertz.cu
	$(NVCC) $(NVCCFLAGS) -D GPU_TIMER \
    op_par.cu ../common/cutil.cpp \
	  ../common/cmd_arg_reader.cpp ../common/stopwatch.cpp ../common/bank_checker.cpp ../common/stopwatch_linux.cpp \
    $(HERTZ_COMMONDIR)/libcommon.so \
    rose_hertz.cu -o op2_gpu_timer
	$(NVCC) $(NVCCFLAGS) -lrt $^ -o op2_posix_timer

src2src: clean
	@if [ -z "${HERTZ_COMMONDIR}" ]; then echo "ERROR: set HERTZ_COMMONDIR env variable first."; exit 1; fi
	$(OP2ROSE) $(INCLUDEFLAGS) hertz.cpp

clean: 
	rm -f op_par.cu kernels.h res_kernel.cu update_kernel.cu rose_hertz.cpp \
	op2_cpu_timer op2_gpu_timer op2_posix_timer rose_hertz.cu blank.cpp
