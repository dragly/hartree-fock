TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += src

!nogui {
    SUBDIRS += hfgui
}
!notools {
    SUBDIRS += tools
}
!notests {
    SUBDIRS += tests
}

BUILD_DIR = $$OUT_PWD
