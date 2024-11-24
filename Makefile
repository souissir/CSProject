# Compiler and linker settings
CC = g++  # Compiler C++
CFLAGS = -g -Wall -Wextra -std=c++14 -Isrc  # Inclure le répertoire src
LDFLAGS = -larmadillo -llapack  # Lier Armadillo et LAPACK (ajoutez d'autres bibliothèques si nécessaire)

# Directories
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

# Targets
TARGET = $(BIN_DIR)/main  # L'exécutable principal (main.cpp + solver.cpp)

# Source and object files
SRC_FILES = $(SRC_DIR)/main.cpp $(SRC_DIR)/solver.cpp $(SRC_DIR)/greeks.cpp # Inclure main.cpp et solver.cpp
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))

# Create the build and bin directories if not present
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)  # Création explicite du dossier bin

# Compile object files for the main application
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Link the main executable
$(TARGET): $(OBJ_FILES) | $(BIN_DIR)  # Ajouter une dépendance pour s'assurer que $(BIN_DIR) existe
	$(CC) $(CFLAGS) $(OBJ_FILES) -o $@ $(LDFLAGS)

# Rule to compile the program
.PHONY: make
make: $(TARGET)

# Rule to run the program
.PHONY: run
run: $(TARGET)
	./$(TARGET)

# Clean build and bin directories
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) $(TARGET)
