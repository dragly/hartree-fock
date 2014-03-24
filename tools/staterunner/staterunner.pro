include(../../defaults.pri)
include(../../defaults-nonsrc.pri)

TEMPLATE = app
CONFIG -= app_bundle
CONFIG -= qt

SOURCES = main.cpp

LIBS += -lhdf5 -lhdf5_cpp

OTHER_FILES += \
    testconfig.cfg
