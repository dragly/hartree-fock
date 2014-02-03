TEMPLATE = subdirs
CONFIG += ordered gui
SUBDIRS += libs app
gui {
    SUBDIRS += hfgui
}
