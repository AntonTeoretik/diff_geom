TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        algebra.cpp \
        main.cpp \
        manifold.cpp \
        metric.cpp

HEADERS += \
    algebra.h \
    constants.h \
    manifold.h \
    metric.h
