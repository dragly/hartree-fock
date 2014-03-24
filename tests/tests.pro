include(../defaults.pri)
include(../defaults-nonsrc.pri)

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
    gaussiandensity.cpp \
    parser.cpp \
    systems.cpp

LIBS += -lunittest++ -L../src/libs -lhartree-fock

HEADERS += \
    basisfunction.h
