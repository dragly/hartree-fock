TEMPLATE = subdirs
CONFIG += ordered
SUBDIRS += src

!noapp {
    SUBDIRS += app
    app.depends = src
}
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

BUILD_DIR = $$OUT_PWD
