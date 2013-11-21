TEMPLATE = subdirs
CONFIG += ordered
SUBDIRS += libs app
gui {
    SUBDIRS += hfgui
}
