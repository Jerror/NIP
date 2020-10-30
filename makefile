SHELL=zsh

WARN = -D"SUPPRESS_NIP_WARNINGS"

INC_DIR = ./
DEBUG = -Og -g3 -fno-omit-frame-pointer -U"NDEBUG"
NDEBUG = -Ofast -funroll-loops -D"NDEBUG"
CXXFLAGS_ADD = -march=native
CXXFLAGS = -Wall -I$(INC_DIR) $(CXXFLAGS_ADD) $(NDEBUG) $(WARN) -pipe
#-Winline 

all: main

main: main.cpp NIP.hpp
	g++ $(CXXFLAGS) -o $@ main.cpp 
