include(../defaults.pri)
include(../defaults-nonsrc.pri)

TARGET = hartree-fock

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -lyaml-cpp

OTHER_FILES += \
    configs/H2O.yaml \
    configs/H2.yaml \
    configs/O2.yaml \
    configs/CH4.yaml
