TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        algebra.cpp \
        color.cpp \
        ellipsoid.cpp \
        genpoints.cpp \
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
    genpoints.h \
    graphics.h \
    manifold.h \
    metric.h \
    weights.h

QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS += -O3

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

