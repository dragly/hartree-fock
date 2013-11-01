include(../../defaults.pri)

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
    math/gaussiantypeoverlapintegral.cpp \
    math/gaussiantypekineticintegral.cpp \
    math/gaussiantypecoloumbattractionintegral.cpp \
    math/gaussiantypeelectroninteractionintegral.cpp

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
    math/gaussiantypeoverlapintegral.h \
    math/gaussiantypekineticintegral.h \
    math/gaussiantypecoloumbattractionintegral.h \
    math/gaussiantypeelectroninteractionintegral.h

OTHER_FILES += \
    defaults.pri
