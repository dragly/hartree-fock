LIBS += -larmadillo -llapack -lblas
LIBS += -lboost_regex

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS

CURRENT_COMPILER = $$QMAKE_CXX
QMAKE_CXX = ccache $$CURRENT_COMPILER

INCLUDEPATH += $$TOP_PWD/src/libs
SRC_DIR = $$TOP_PWD

copydata.commands = $(COPY_DIR) $$TOP_PWD/utils/* $$OUT_PWD
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata


folder_01.source = $$TOP_PWD/utils
folder_01.target = utils

DEPLOYMENTFOLDERS = folder_01

for(deploymentfolder, DEPLOYMENTFOLDERS) {
    item = item$${deploymentfolder}
    greaterThan(QT_MAJOR_VERSION, 4) {
        itemsources = $${item}.files
    } else {
        itemsources = $${item}.sources
    }
    $$itemsources = $$eval($${deploymentfolder}.source)
    itempath = $${item}.path
    $$itempath= $$eval($${deploymentfolder}.target)
    export($$itemsources)
    export($$itempath)
    DEPLOYMENT += $$item
}
