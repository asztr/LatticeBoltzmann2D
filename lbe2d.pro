 QT          += opengl

 HEADERS     = glwidget.h \
               window.h \
               Mouse.h \
               Lattice.h \
               Cell.h \
               Common.h \
               fftw++/Array.h \
               fftw++/fftw++.h

 SOURCES     = glwidget.cpp \
               main.cpp \
               window.cpp \
               Mouse.cpp \
               Lattice.cpp \
               Cell.cpp \
               fftw++/fftw++.cc

QMAKE_CXXFLAGS += -std=c++11 -Wno-narrowing

#optimizaciones y openMP
#QMAKE_CXXFLAGS += -03 -fopenmp -funroll-loops -ffast-math -fomit-frame-pointer
#QMAKE_LFLAGS += -fopenmp

#profiler
#QMAKE_CXXFLAG += -std=c++0x -pg -g
#QMAKE_LFLAGS += -pg

LIBS += -lfftw3


 # install
 target.path = $$[QT_INSTALL_EXAMPLES]/opengl/2dpainting
 sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS lbe2d.pro
 sources.path = $$[QT_INSTALL_EXAMPLES]/opengl/2dpainting
 INSTALLS += target sources

 symbian: include($$QT_SOURCE_TREE/examples/symbianpkgrules.pri)
