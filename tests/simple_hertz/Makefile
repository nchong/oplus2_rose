COMMONDIR := ../common
OP2ROSE := ../../translator/op2rose

NVCC := /data/cuda30/cuda/bin/nvcc
INCLUDEFLAGS := -I$(COMMONDIR)
NVCCFLAGS := -arch sm_13 -g $(INCLUDEFLAGS) -D OP_DIAGS=0 -deviceemu

VPATH = $(COMMONDIR)

SOURCES = rose_hertz.cpp op_par.cu cutil.cpp \
	cmd_arg_reader.cpp stopwatch.cpp bank_checker.cpp stopwatch_linux.cpp

all: src2src build

build: $(SOURCES)
	$(NVCC) $(NVCCFLAGS) $^ -o simple_hertz

src2src: clean
	$(OP2ROSE) $(INCLUDEFLAGS) hertz.cpp

clean: 
	rm -f op_par.cu kernels.h res_kernel.cu rose_hertz.cpp \
	simple_hertz rose_hertz.cu blank.cpp
