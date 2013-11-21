#include <QtGui/QGuiApplication>
#include <hartreefock.h>
#include <volumeeffect.h>
#include "qtquick2applicationviewer.h"

int main(int argc, char *argv[])
{
    qmlRegisterType<HartreeFock>("Dragly", 1, 0, "HartreeFock");
    qmlRegisterType<VolumeShaderProgram>("Dragly", 1, 0, "VolumeShaderProgram");
    QGuiApplication app(argc, argv);

    QtQuick2ApplicationViewer viewer;
    viewer.setMainQmlFile(QStringLiteral("qml/hfgui/main.qml"));
    viewer.showExpanded();

    return app.exec();
}
