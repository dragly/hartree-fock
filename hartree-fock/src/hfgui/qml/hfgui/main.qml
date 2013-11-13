import QtQuick 2.0
import Qt3D 2.0
import Qt3D.Shapes 2.0
import Dragly 1.0
import QtQuick.Controls 1.0
import QtQuick.Dialogs 1.0
import QtQuick.Layouts 1.0

Rectangle {
    width: 1280
    height: 800

    RowLayout {
        anchors.fill: parent
        anchors.margins: 5
        Viewport {
            id: mainViewport
            property alias multiplier: mainDensityPlotter.multiplier
            property alias useSquareRootDensity: mainDensityPlotter.useSquareRootDensity
            property alias volumeShaderQuality: mainDensityPlotter.quality
            fillColor: "black"

            Layout.fillWidth: true
            Layout.fillHeight: true

            HartreeFock {
                id: hartreeFock
            }

            Item3D {
                id: volumeItem
                property double maxMinDifference: hartreeFock.voxelEdgeMax - hartreeFock.voxelEdgeMin
                cullFaces: Item3D.CullBackFaces
                scale: maxMinDifference
                mesh: Mesh {
                    source: "cube.obj"
                }
                effect: VolumeShaderProgram {
                    id: mainDensityPlotter
                    property real multiplier: densityMultiplierSlider.value
                    property bool useSquareRootDensity: useSquareRootDensityCheckBox.checked
                    property real quality: volumeShaderQualitySlider.value
                    blending: true
                    vertexShaderSource: "scalarvolume.vert"
                    fragmentShaderSource: "scalarvolume.frag"
                    positionReader: hartreeFock
                }
            }
        }
        ColumnLayout {
            Label {
                text: qsTr("Density multiplier:")
            }
            Slider {
                id: densityMultiplierSlider
                Layout.minimumWidth: 200
                minimumValue: 1e-6
                maximumValue: 100.0
                value: 30
            }
            CheckBox {
                id: useSquareRootDensityCheckBox
                checked: false
                text: "Use sqrt(densityValue)"
            }
            Label {
                text: qsTr("Volume shader quality:")
            }
            Slider {
                id: volumeShaderQualitySlider
                Layout.minimumWidth: 200
                minimumValue: 10
                maximumValue: 1000
                value: 100
            }
            Label {
                text: qsTr("Energy:") + " " + hartreeFock.energy.toFixed(2)
            }

            Item {
                Layout.fillHeight: true
            }
        }
    }


}
