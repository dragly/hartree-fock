include(../defaults.pri)
include(../defaults-nonsrc.pri)

QT += 3dquick

# Add more folders to ship with the application, here
folder_01.source = qml/hfgui
folder_01.target = qml
DEPLOYMENTFOLDERS = folder_01

# Additional import path used to resolve QML modules in Creator's code model
QML_IMPORT_PATH =

# If your application uses the Qt Mobility libraries, uncomment the following
# lines and add the respective components to the MOBILITY variable.
# CONFIG += mobility
# MOBILITY +=

# The .cpp file which was generated for your project. Feel free to hack it.
SOURCES += main.cpp \
    volumeeffect.cpp \
    hartreefock.cpp \
    settings.cpp

# Installation path
# target.path =

# Please do not modify the following two lines. Required for deployment.
include(qtquick2applicationviewer/qtquick2applicationviewer.pri)
qtcAddDeployment()

HEADERS += \
    volumeeffect.h \
    volumeeffect_p.h \
    hartreefock.h \
    settings.h

OTHER_FILES += \
    qml/hfgui/scalarvolume.frag \
    defaults-nonsrc.pri

LIBS += -L$$PWD/src -lhartree-fock
LIBS += -lyaml-cpp

#copydata.commands = $(COPY_DIR) ../data/* $$OUT_PWD
first.depends = $(first) copydeploymentfolders copydata
export(first.depends)
#export(copydata.commands)
#QMAKE_EXTRA_TARGETS += first copydata copydeploymentfolders
