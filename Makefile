EIGEN := $(shell pkg-config eigen3 --cflags)
ROOT := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)

SRC := $(wildcard src/*.cpp)

C_FLAGS := -std=c++17
C_FLAGS_ADD := -Wall -Werror -Wextra -O2

#hw3: proj3.exe proj2.exe proj4.exe
    

%.exe: main/%.cpp
	@echo "COMPILING $< --> $@"
	g++ $(C_FLAGS) $< -o bin/$@ -I inc $(SRC) $(ROOT) $(ROOT_LIBS)

dump:
	@echo $(EIGEN)
	@echo $(ROOT)
	@echo $(ROOT_LIBS)
	@echo $(SRC)
