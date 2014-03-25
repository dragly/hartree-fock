QMAKE_CXX = mpicxx

include(../../defaults.pri)
include(../../defaults-nonsrc.pri)

TEMPLATE = app
CONFIG -= app_bundle
CONFIG -= qt

SOURCES = main.cpp

LIBS += -lhdf5 -lhdf5_cpp
LIBS += -lboost_filesystem -lboost_system -lboost_mpi -lboost_serialization
LIBS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)

OTHER_FILES += \
    testconfig.cfg
