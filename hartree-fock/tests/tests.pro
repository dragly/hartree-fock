include(../defaults.pri)

TARGET = hartree-fock-tests

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    hydrogenplot.cpp \
    helium-tests.cpp

LIBS += -lunittest++ -L$$TOP_OUT_PWD/src/libs -lhartree-fock
