TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        algebra.cpp \
        atlas.cpp \
        geometry.cpp \
        main.cpp

HEADERS += \
    algebra.h \
    atlas.h \
    constants.h \
    geometry.h
