include($$TOP_PWD/defaults.pri)

TARGET = hartree-fock-tests

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    helium-tests.cpp \
    gaussiantypeintegrals.cpp \
    boysfunction.cpp \
    hydrogen.cpp \
    gaussianoxygen.cpp \
    gaussiandensity.cpp

LIBS += -lunittest++ -L$$TOP_OUT_PWD/src/libs -lhartree-fock

HEADERS += \
    basisfunction.h
