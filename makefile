CXX = g++
CXXFLAGS = -std=c++17 -Wshadow -Wextra -lglfw -lGLEW -lGLU -lGL -larmadillo -llapack -lblas

BIN = bin
SRC = src
INCLUDE = include
LIB = lib
EXECUTABLE = SNAC


all: $(BIN)/$(EXECUTABLE)

run: clean all
	clear
	@echo "Executing."
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	@echo "Building."
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@

clean:
	@echo "Clearing."
	-rm $(BIN)/*
