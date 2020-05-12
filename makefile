CXX = g++
CC  = $(CXX)

CXXFLAGS = -std=c++1z

# Partie commentÃ©e : choisissez les options que vous voulez avoir
#                    en dÃ©commentant la/les lignes correspondantes
#
CXXFLAGS += -pedantic -Wall -Wextra # pour les purs et durs
CXXFLAGS += -g #                    # pour debugger
# CXXFLAGS += -pg #                 # pour profiler
# LDFLAGS  += -pg #                 # pour profiler
# CXXFLAGS += -O2 #                 # pour optimiser la vitesse

LIB=LIB/
SRC=SRC/
BIN=

all: $(BIN)PIMC $(BIN)PIMC_H3part


$(LIB)PIMC.o : $(SRC)PIMC.cc
	$(CC) $(CXXFLAGS) -c -o $@ $<

$(LIB)PIMC_H3part.o : $(SRC)PIMC_H3part.cc
	$(CC) $(CXXFLAGS) -c -o $@ $<



$(BIN)PIMC: $(LIB)PIMC.o
	$(CC) $(CXXFLAGS) $(LIB)PIMC.o -o $(BIN)PIMC


$(BIN)PIMC_H3part: $(LIB)PIMC_H3part.o
	$(CC) $(CXXFLAGS) $(LIB)PIMC_H3part.o -o $(BIN)PIMC_H3part
