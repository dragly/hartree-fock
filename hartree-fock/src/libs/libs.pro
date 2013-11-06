include($$TOP_PWD/defaults.pri)

TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = hartree-fock

SOURCES += hartreesolver.cpp \
    math/vector3.cpp \
    electronsystems/helium/heliumhartree.cpp \
    electronsystems/hydrogen/hydrogenmolecule.cpp \
    hartreefocksolver.cpp \
    electronsystems/hydrogen/multihydrogen.cpp \
    electronsystems/electronsystem.cpp \
    basisfunctions/gaussiantypeorbitalintegrator.cpp \
    basisfunctions/gaussiantypeorbital.cpp \
    basisfunctions/atomicbasisfunction.cpp \
    math/boysfunction.cpp \
    math/boysfunctionintermediate.cpp \
    hermiteintegral.cpp \
    math/hermiteexpansioncoefficient.cpp \
    basisfunctions/gaussiantypeorbital/integrals/gaussiantypeelectroninteractionintegral.cpp \
    basisfunctions/gaussiantypeorbital/integrals/gaussiantypekineticintegral.cpp \
    basisfunctions/gaussiantypeorbital/integrals/gaussiantypeoverlapintegral.cpp \
    basisfunctions/gaussiantypeorbital/integrals/gaussiantypecoloumbattractionintegral.cpp

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
    basisfunctions/atomicbasisfunction.h\
    math/boysfunction.h \
    math/boysfunctionintermediate.h \
    hermiteintegral.h \
    math/hermiteexpansioncoefficient.h \
    basisfunctions/gaussiantypeorbital/integrals/gaussiantypecoloumbattractionintegral.h \
    basisfunctions/gaussiantypeorbital/integrals/gaussiantypeelectroninteractionintegral.h \
    basisfunctions/gaussiantypeorbital/integrals/gaussiantypekineticintegral.h \
    basisfunctions/gaussiantypeorbital/integrals/gaussiantypeoverlapintegral.h

OTHER_FILES += \
    defaults.pri \
    ../../.qmake.conf
