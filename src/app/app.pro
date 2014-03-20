include($$TOP_PWD/defaults.pri)

TEMPLATE = app
CONFIG -= app_bundle
CONFIG -= qt

SOURCES = main.cpp

LIBS += -lhartree-fock -L../libs
LIBS += -lconfig++
LIBS += -lhdf5 -lhdf5_cpp

OTHER_FILES += \
    testconfig.cfg
