#INCDIR = -I$(ROOTSYS)/include -Iinc
INCDIR = -Iinc -I$(LT)/include

LIBS = -lstdc++ -L/usr/local/lib -lm -L$(LT)/lib64 -looptools -lgfortran -lm -lquadmath -lm #$(shell root-config --evelibs) -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO -lMLP -lm -lRooFit -lRooFitCore
CC = fcc -fPIC #g++ -v

CFLAGS = -g -ansi -pedantic -W -Wall -Wextra -std=c++11 -m64 -fPIC #-std=c++17
LDFLAGS= -shared
#LDFLAGS= -L$(LT)/lib64 -looptools

TARGET = libloopintegral.so

SRC = $(wildcard src/*.cxx)

OBJECTS = $(subst src,obj,$(patsubst %.cxx,%.o,$(SRC)))

HEADERS = $(wildcard $(INCDIR)/*.h)

CONDA_PATH = /home/andres/anaconda3/envs/higgs
CONDA_ACTIVATE = $(shell source activate $(CONDA_PATH))
CONDA_ENV = higgs

.PHONY: all clean

all : $(TARGET)
	@echo "DONE!"

obj/%.o: src/%.cxx $(HEADERS)
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(TARGET): $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(LIBS) $(OBJECTS) -o $@ $(LDFLAGS)
obj:
	mkdir -p $@
clean:
	rm -rf ./obj
	rm -rf $(TARGET)
