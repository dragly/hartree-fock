#include <QtGui/QGuiApplication>
#include <hartreefock.h>
#include <volumeeffect.h>
#include "qtquick2applicationviewer.h"
#include "settings.h"

int main(int argc, char *argv[])
{
    QCoreApplication::setOrganizationName("dragly");
    QCoreApplication::setOrganizationDomain("dragly.org");
    QCoreApplication::setApplicationName("Denseness");

    qmlRegisterType<HartreeFock>("Dragly", 1, 0, "HartreeFock");
    qmlRegisterType<VolumeShaderProgram>("Dragly", 1, 0, "VolumeShaderProgram");
    qmlRegisterType<Settings>("Settings", 1, 0, "Settings");

    QGuiApplication app(argc, argv);

    QtQuick2ApplicationViewer viewer;
    viewer.setMainQmlFile(QStringLiteral("qml/hfgui/main.qml"));
    viewer.showExpanded();

    return app.exec();
}
