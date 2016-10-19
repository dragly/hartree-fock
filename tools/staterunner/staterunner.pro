QMAKE_CXX = mpicxx

include(../../defaults.pri)
include(../../defaults-nonsrc.pri)

TEMPLATE = app
CONFIG -= app_bundle
CONFIG -= qt

SOURCES = main.cpp

LIBS += -lhdf5 -lhdf5_cpp
LIBS += -lboost_filesystem -lboost_system -lboost_mpi -lboost_serialization
LIBS += -lyaml-cpp
LIBS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)

INCLUDEPATH += /usr/include/hdf5/serial/
LIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial/

OTHER_FILES += \
    testconfig.cfg
