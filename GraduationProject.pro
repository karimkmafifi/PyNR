#-------------------------------------------------
#
# Project created by QtCreator 2017-05-26T06:04:34
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = GraduationProject
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += "/Users/karimafifi/Desktop/PyNR-master/includelibs"

SOURCES += main.cpp\
        mainwindow.cpp \
    score.cpp \
    parser.cpp \
    ELEMENT_DATA_new_GN.cpp \
    mw_match_GN.cpp \
    prepare_mol.cpp \
    scoring_terms.cpp \
    linearRegression.cpp \
    getxswindow.cpp \
    trainwindow.cpp \
    vswindow.cpp \
    randomforest.cpp \
    scorepdbbind.cpp \
    scorevs.cpp \
    testwindow.cpp \
    testpdbbind.cpp

HEADERS  += mainwindow.h \
    common.h \
    input.h \
    atom.h \
    atom_constants.h \
    filesystem.h \
    gss_sphere_points.h \
    parser.h \
    scoring_terms.h \
    score.h \
    molecule.h \
    bond.h \
    ELEMENT_DATA_new_GN.h \
    mw_match_GN.h \
    octree_GN.hpp \
    stl_ptr_GN.hpp \
    prepare_mol.h \
    linearRegression.h \
    getxswindow.h \
    trainwindow.h \
    vswindow.h \
    randomforest.h \
    scorepdbbind.h \
    scorevs.h \
    testwindow.h \
    testpdbbind.h

FORMS    += mainwindow.ui \
    getxswindow.ui \
    trainwindow.ui \
    vswindow.ui \
    testwindow.ui

RESOURCES += \
    mainwindow.qrc

DISTFILES +=
