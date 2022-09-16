TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        algebra.cpp \
        color.cpp \
        ellipsoid.cpp \
        main.cpp \
        manifold.cpp \
        metric.cpp

HEADERS += \
    algebra.h \
    color.h \
    constants.h \
    ellipsoid.h \
    manifold.h \
    metric.h
