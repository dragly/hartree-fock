TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    hartreesolver.cpp \
    math/vector3.cpp \
    basisfunctions/basisfunction.cpp \
    basisfunctions/helium/heliumhartree.cpp \
    basisfunctions/hydrogen/hydrogenmolecule.cpp \
    hartreefocksolver.cpp \
    basisfunctions/hydrogen/multihydrogen.cpp

HEADERS += \
    hartreesolver.h \
    math/vector3.h \
    basisfunctions/basisfunction.h \
    basisfunctions/helium/heliumhartree.h \
    basisfunctions/hydrogen/hydrogenmolecule.h \
    hartreefocksolver.h \
    basisfunctions/hydrogen/multihydrogen.h

LIBS += -larmadillo -llapack -lblas

include(defaults.pri)

OTHER_FILES += \
    defaults.pri
