#-------------------------------------------------
#
# Project created by QtCreator 2015-09-23T08:32:31
#
#-------------------------------------------------

QT       -= core
QT       -= gui

OBJECTS_DIR = build
DESTDIR = bin

TARGET = openmp_performance_test
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

HEADERS += qute.h
SOURCES += qute.cpp

SOURCES += main.cpp

QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp -std=c++11
LIBS += -fopenmp -lfftw3 -lfftw3_threads
