#-------------------------------------------------
#
# Project created by QtCreator 2018-05-28T12:31:59
#
#-------------------------------------------------

QT          += core
QT          += gui
QT          += opengl
#QT          += openglextensions
QT          += charts

greaterThan(QT_MAJOR_VERSION, 5): QT += widgets

TARGET       = hpg_dhunter
TEMPLATE     = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES     += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES     += main.cpp \
               hpg_dhunter.cpp \
               files_worker.cpp \
               refgen.cpp

HEADERS     += \
               data_pack.h \
               hpg_dhunter.h \
               files_worker.h \
               refgen.h

FORMS       += \
               hpg_dhunter.ui

OTHER_FILES +=

CONFIG      += C++11

DESTDIR      = $$system(pwd)
OBJECTS_DIR  = $$DESTDIR/Obj

#----------------------------------------------------------------------
#-----------------------------cuda settings----------------------------
#----------------------------------------------------------------------
# cuda sources
CUDA_SOURCES += haar_v10.cu

# path to cuda sdk installation
CUDA_DIR      = /usr/local/cuda

# path to header and libs files
INCLUDEPATH  += $$CUDA_DIR/include
QMAKE_LIBDIR += $$CUDA_DIR/lib64

# cuda architecture
#CUDA_ARCH    = sm_35       # minimum compute capability (version) for dynamic parallelism feature support
CUDA_ARCH     = sm_61
#CUDA_ARCH    = sm_50
#CUDA_ARCH    = sm_52

# libs used in the code
LIBS         += -lcudart -lcuda -lcudadevrt

# some nvcc compiler flags
NVCCFLAGS     = --compiler-options \
                -fno-strict-aliasing \
                -std=c++11 \
                -use_fast_math \
                --ptxas-options=-v

# prepare the extra compiler configuration
CUDA_INC      = $$join(INCLUDEPATH,' -I','-I',' ')

# prepare intermediate CUDA compiler  - - - - - - - - - - - - - - - - - - - - - - - -
# this is neccesary because there are more than one __global__ functions
# then it must be compiled with dynamic parallelism
cudaIntr.input  = CUDA_SOURCES
cudaIntr.output = ${OBJECTS_DIR}${QMAKE_FILE_BASE}.o

# tweak arch according to your hws compute capability
cudaIntr.commands = $$CUDA_DIR/bin/nvcc \
                    -m64 \                  # type of machine
                    -g \                    # debug mode for host code
                    -G \                    # debug mode for device code
                    -arch=$$CUDA_ARCH \     # device architecture for files of data
                    -dc \                   # dynamic parallelism compiler
                    $$NVCCFLAGS \
                    $$CUDA_INC \
                    $$LIBS \
                    ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT}

# set our variable out.
# these obj files need to be used to create the link obj file
# and used in our final gcc compilation
cudaIntr.variable_out  = CUDA_OBJ
cudaIntr.variable_out += OBJECTS
cudaIntr.clean         = cudaIntrObj/*.o

# tell Qt that we want add more stuff to the Makefile
QMAKE_EXTRA_COMPILERS += cudaIntr

# prepare the linking compiler step - - - - - - - - - - - - - - - - - - - - - - - - -
cuda.input    = CUDA_OBJ
cuda.output   = ${QMAKE_FILE_BASE}_link.o

# Tweak arch according to your hws compute capability
cuda.commands        = $$CUDA_DIR/bin/nvcc \
                       -m64 \
                       -g \
                       -G \
                       -arch=$$CUDA_ARCH \
                       -dlink ${QMAKE_FILE_NAME} \
                       -o ${QMAKE_FILE_OUT}

cuda.dependency_type = TYPE_C

cuda.depend_command  = $$CUDA_DIR/bin/nvcc \
                       -g \
                       -G \
                       -M \                     # link the previous object with main exec
                       $$CUDA_INC \
                       $$NVCCFLAGS \
                       ${QMAKE_FILE_NAME}

# tell Qt that we want add more stuff to the Makefile
QMAKE_EXTRA_COMPILERS += cuda

DISTFILES += \
    haar_v10.cu

RESOURCES += \
    recursos.qrc
