TEMPLATE = subdirs

!nogui {
    SUBDIRS += hfgui
    hfgui.depends = src
}
!notools {
    SUBDIRS += tools
    tools.depends = src
}
!notests {
    SUBDIRS += tests
    tests.depends = src
}

SUBDIRS += src

BUILD_DIR = $$OUT_PWD
