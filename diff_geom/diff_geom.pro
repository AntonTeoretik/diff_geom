TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        algebra.cpp \
        color.cpp \
        ellipsoid.cpp \
        graphics.cpp \
        main.cpp \
        manifold.cpp \
        metric.cpp

HEADERS += \
    algebra.h \
    bitimage/bitmap_image.hpp \
    color.h \
    constants.h \
    ellipsoid.h \
    graphics.h \
    manifold.h \
    metric.h

QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS += -O3
