include(../defaults.pri)

TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = hartree-fock

SOURCES += \
    math/vector3.cpp \
    electronsystems/helium/heliumhartree.cpp \
    electronsystems/hydrogen/hydrogenmolecule.cpp \
    electronsystems/hydrogen/multihydrogen.cpp \
    electronsystems/electronsystem.cpp \
    basisfunctions/gaussiantypeorbital.cpp \
    basisfunctions/atomicbasisfunction.cpp \
    math/boysfunction.cpp \
    math/boysfunctionintermediate.cpp \
    hermiteintegral.cpp \
    math/hermiteexpansioncoefficient.cpp \
    basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.cpp \
    basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.cpp \
    basisfunctions/gaussian/integrals/gaussiankineticintegral.cpp \
    basisfunctions/gaussian/integrals/gaussianoverlapintegral.cpp \
    basisfunctions/gaussian/gaussiancontractedorbital.cpp \
    basisfunctions/gaussian/gaussianprimitiveorbital.cpp \
    electronsystems/gaussian/gaussiansystem.cpp \
#    electronsystems/gaussian/gaussiannitrogen431g.cpp \
#    electronsystems/gaussian/gaussianoxygen431g.cpp \
    parsers/turbomoleparser.cpp \
    electronsystems/gaussian/gaussiancore.cpp \
    hf.cpp \
    solvers/unrestrictedhartreefocksolver.cpp \
    solvers/hartreefocksolver.cpp \
    solvers/hartreesolver.cpp \
    solvers/restrictedhartreefocksolver.cpp

HEADERS += \
    math/vector3.h \
    electronsystems/electronsystem.h \
    electronsystems/helium/heliumhartree.h \
    electronsystems/hydrogen/hydrogenmolecule.h \
    electronsystems/hydrogen/multihydrogen.h \
    basisfunctions/gaussiantypeorbital.h \
    basisfunctions/atomicbasisfunction.h\
    math/boysfunction.h \
    math/boysfunctionintermediate.h \
    hermiteintegral.h \
    math/hermiteexpansioncoefficient.h \
    basisfunctions/gaussian/integrals/gaussiancoloumbattractionintegral.h \
    basisfunctions/gaussian/integrals/gaussianelectroninteractionintegral.h \
    basisfunctions/gaussian/integrals/gaussiankineticintegral.h \
    basisfunctions/gaussian/integrals/gaussianoverlapintegral.h \
    basisfunctions/gaussian/gaussiancontractedorbital.h \
    basisfunctions/gaussian/gaussianprimitiveorbital.h \
    electronsystems/gaussian/gaussiansystem.h \
#    electronsystems/gaussian/gaussiannitrogen431g.h \
#    electronsystems/gaussian/gaussianoxygen431g.h \
    parsers/turbomoleparser.h \
    electronsystems/gaussian/gaussiancore.h \
    hf.h \
    solvers/unrestrictedhartreefocksolver.h \
    solvers/restrictedhartreefocksolver.h \
    solvers/hartreefocksolver.h \
    solvers/hartreesolver.h

OTHER_FILES += \
    defaults.pri \
    ../../.qmake.conf
