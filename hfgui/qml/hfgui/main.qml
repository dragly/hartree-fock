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
//            property alias useSquareRootDensity: mainDensityPlotter.useSquareRootDensity
            property alias volumeShaderQuality: mainDensityPlotter.quality
            fillColor: backgroundColorPicker.color
            fovzoom: false
            camera: Camera {
                center: Qt.vector3d(0,0,0)
                eye: Qt.vector3d(10, 0, 0)
                nearPlane: 0.5
                farPlane: 100
                fieldOfView: 60
            }

            Layout.fillWidth: true
            Layout.fillHeight: true

            HartreeFock {
                id: hartreeFock
                orbital: allOrbitalsCheckBox.checked ? -1 : orbitalSlider.value - 1
            }

            Item3D {
                id: volumeItem
                property double maxMinDifference: hartreeFock.voxelEdgeMax - hartreeFock.voxelEdgeMin
                cullFaces: Item3D.CullFrontFaces
                scale: maxMinDifference
                mesh: Mesh {
                    source: "cube.obj"
                }
                effect: VolumeShaderProgram {
                    id: mainDensityPlotter
                    property real multiplier: densityMultiplierSlider.expValue
//                    property bool useSquareRootDensity: useSquareRootDensityCheckBox.checked
                    property real contrast: densityContrastSlider.expValue
                    property real quality: volumeShaderQualitySlider.value
                    property color standardColor: standardColorPicker.color
                    property color highlightColor: highlightColorPicker.color
//                    onEffectChanged:
                    blending: true
                    vertexShaderSource: "scalarvolume.vert"
                    fragmentShaderSource: "scalarvolume.frag"
                    positionReader: hartreeFock
                }
            }
        }
        ColumnLayout {
            Layout.minimumWidth: 200
            spacing: 5
            Label {
                text: qsTr("Orbital:") + " " + orbitalSlider.value
            }
            Slider {
                id: orbitalSlider
                Layout.preferredWidth: parent.width
                minimumValue: 1
                maximumValue: hartreeFock.orbitalCount
                stepSize: 1
                enabled: !allOrbitalsCheckBox.checked
                opacity: allOrbitalsCheckBox.checked ? 0.3 : 1.0
            }
            CheckBox {
                id: allOrbitalsCheckBox
                text: "All"
                checked: true
            }

            Label {
                text: qsTr("Density multiplier:")
            }
            Slider {
                id: densityMultiplierSlider
                property real expValue: Math.exp(value * densityContrastSlider.expValue*maximumValue) // Makes the slider scale logarithmic
                Layout.preferredWidth: parent.width
                minimumValue: -3
                maximumValue: 3
                value: 2
            }
            Label {
                text: qsTr("Density contrast:")
            }
            Slider {
                id: densityContrastSlider
                property real expValue: Math.exp(value)
                Layout.preferredWidth: parent.width
                minimumValue: -3
                maximumValue: 3
                value: 0
            }
//            CheckBox {
//                id: useSquareRootDensityCheckBox
//                checked: false
//                text: "Use sqrt(densityValue)"

//            }
            Label {
                text: qsTr("Volume shader quality:")
            }
            Slider {
                id: volumeShaderQualitySlider
                Layout.preferredWidth: parent.width
                minimumValue: 1e1
                maximumValue: 1e3
                value: 3e2
            }
            Label {
                text: qsTr("Energy:") + " " + hartreeFock.energy.toFixed(2)
            }
            Label {
                text: "Color:"
            }
            ColorPicker {
                id: standardColorPicker
                color: "#08306B"
                Layout.preferredWidth: parent.width
                height: 30
            }

            Label {
                text: "Highlight:"
            }
            ColorPicker {
                id: highlightColorPicker
                color: "#44C5DAEE"
                Layout.preferredWidth: parent.width
                height: 30
            }

            Label {
                text: "Background:"
            }
            ColorPicker {
                id: backgroundColorPicker
                color: "white"
                Layout.preferredWidth: parent.width
                height: 30
            }

            Item {
                Layout.fillHeight: true
            }
        }
    }


}
