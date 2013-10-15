TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    hartreesolver.cpp \
    math/vector3.cpp \
    electronsystems/helium/heliumhartree.cpp \
    electronsystems/hydrogen/hydrogenmolecule.cpp \
    hartreefocksolver.cpp \
    electronsystems/hydrogen/multihydrogen.cpp \
    electronsystems/electronsystem.cpp \
    basisfunctions/gaussiantypeorbitalintegrator.cpp \
    basisfunctions/gaussiantypeorbital.cpp \
    basisfunctions/atomicbasisfunction.cpp

HEADERS += \
    hartreesolver.h \
    math/vector3.h \
    electronsystems/electronsystem.h \
    electronsystems/helium/heliumhartree.h \
    electronsystems/hydrogen/hydrogenmolecule.h \
    hartreefocksolver.h \
    electronsystems/hydrogen/multihydrogen.h \
    basisfunctions/gaussiantypeorbitalintegrator.h \
    basisfunctions/gaussiantypeorbital.h \
    basisfunctions/atomicbasisfunction.h

LIBS += -larmadillo -llapack -lblas

include(defaults.pri)

OTHER_FILES += \
    defaults.pri
