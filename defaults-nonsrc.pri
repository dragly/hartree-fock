TOP_OUT_PWD = $$shadowed($$PWD)
LIBS += -L$$TOP_OUT_PWD/src -lhartree-fock

INCLUDEPATH += $$PWD/src
SRC_DIR = $$PWD

copydata.commands = $(COPY_DIR) $$PWD/data/* $$OUT_PWD
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata

folder_01.source = $$PWD/data
folder_01.target = data

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
