LIBS += -larmadillo -llapack -lblas

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS

CURRENT_COMPILER = $$QMAKE_CXX
QMAKE_CXX = ccache $$CURRENT_COMPILER

PRO_PATH=$$PWD
INCLUDEPATH += $$PRO_PATH/..
SRC_DIR = $$PRO_PATH
