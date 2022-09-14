TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        algebra.cpp \
        geometry.cpp \
        main.cpp \
        manifold.cpp

HEADERS += \
    algebra.h \
    constants.h \
    geometry.h \
    manifold.h
