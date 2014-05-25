import QtQuick 2.0
import QtQuick.Layouts 1.1
import QtQuick.Dialogs 1.1
import QtQuick.Controls 1.1

Item {
    id: rectRoot
    property color color: "blue"
    property bool _internalChange: false
    Rectangle {
        anchors.fill: parent
        color: rectRoot.color
        opacity: colorDialog.currentAlpha
    }
    Rectangle {
        anchors.fill: parent
        color: "#00000000"
        border.width: 1
        border.color: "#AA000000"
    }
    onColorChanged: {
        if(!_internalChange) {
            colorDialog.currentColor = color
        }
    }
    ColorDialog {
        id: colorDialog
        color: rectRoot.color
        showAlphaChannel: true

        onCurrentColorChanged: {
            _internalChange = true
            rectRoot.color = currentColor
            _internalChange = false
        }

        onAccepted: {
            color = currentColor
        }
        onRejected: {
            currentColor = color
        }
    }

    MouseArea {
        anchors.fill: parent
        onClicked: colorDialog.open()
    }
}
