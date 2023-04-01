QT += widgets

HEADERS += \
        mainwindow.h \
        nodeandedge.h

SOURCES += \
        main.cpp \
        mainwindow.cpp \
        nodeandedge.cpp

RESOURCES = designer.qrc

# install
target.path = $$[QT_INSTALL_EXAMPLES]/widgets/graphicsview/elasticnodes
INSTALLS += target

FORMS += \
    mainwindow.ui
