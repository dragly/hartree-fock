TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    hartreesolver.cpp \
    math/vector3.cpp

HEADERS += \
    hartreesolver.h \
    math/vector3.h

LIBS += -larmadillo -llapack -lblas

include(defaults.pri)

OTHER_FILES += \
    defaults.pri
