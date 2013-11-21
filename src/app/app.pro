include($$TOP_PWD/defaults.pri)

TEMPLATE = app
CONFIG -= app_bundle
CONFIG -= qt

SOURCES = main.cpp

LIBS += -lhartree-fock -L../libs
