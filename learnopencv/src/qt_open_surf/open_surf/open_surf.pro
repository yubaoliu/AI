#-------------------------------------------------
#
# Project created by QtCreator 2012-08-17T21:54:15
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = open_surf
TEMPLATE = app


SOURCES += main.cpp\
        opensurf.cpp \
    ipoint.cpp \
    fasthessian.cpp \
    surf.cpp \
    integral.cpp \
    utils.cpp \
    kmeans.cpp

HEADERS  += opensurf.h \
    ipoint.h \
    kmeans.h \
    utils.h \
    integral.h \
    fasthessian.h \
    surf.h \
    responselayer.h \
    surflib.h

FORMS    += opensurf.ui


INCLUDEPATH +=  C:\Qt\opencv2.4.2\build\include \
                C:\Qt\opencv2.4.2\build\include\opencv \
                C:\Qt\opencv2.4.2\build\include\opencv2

LIBS += C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_core242d.lib    \
        C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_highgui242d.lib  \
        C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_imgproc242d.lib   \
        C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_legacy242.lib   \
        C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_features2d242d.lib  \
        C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_contrib242d.lib     \
        C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_flann242.lib        \
        C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_objdetect242d.lib   \
        C:\Qt\opencv2.4.2\build\x86\vc10\lib\opencv_calib3d242d.lib
