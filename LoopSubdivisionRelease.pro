######################################################################
# Automatically generated by qmake (3.1) Thu Dec 15 03:03:35 2022
######################################################################

QT+=opengl
TEMPLATE = app
TARGET = LoopSubdivisionRelease
INCLUDEPATH += .
LIBS += -lopengl32

# You can make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# Please consult the documentation of the deprecated API in order to know
# how to port your code away from it.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input
HEADERS += ArcBall.h \
           ArcBallWidget.h \
           Cartesian3.h \
           DirectedEdgeSurface.h \
           Homogeneous4.h \
           Matrix4.h \
           Quaternion.h \
           RenderController.h \
           RenderParameters.h \
           RenderWidget.h \
           RenderWindow.h \
           RGBAImage.h \
           RGBAValue.h \
           SphereVertices.h
SOURCES += ArcBall.cpp \
           ArcBallWidget.cpp \
           Cartesian3.cpp \
           DirectedEdgeSurface.cpp \
           Homogeneous4.cpp \
           main.cpp \
           Matrix4.cpp \
           Quaternion.cpp \
           RenderController.cpp \
           RenderWidget.cpp \
           RenderWindow.cpp \
           RGBAImage.cpp \
           RGBAValue.cpp \
           SphereVertices.cpp
