TEMPLATE = subdirs
CONFIG += ordered
SUBDIRS += libs app
!nogui {
    SUBDIRS += hfgui
}
