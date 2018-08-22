CC = g++ -std=c++11
OPTIM_FLAG = -O3
CC_FLAGS = $(OPTIM_FLAG) -Wall -c
LD_FLAGS = $(OPTIM_FLAG) -Wall
LD_LIBS = -lgsl -lgslcblas -lm -pthread
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

all: obj mlev-hiv-sim data

clean:
	rm -f obj/*.o obj/*.d mlev-hiv-sim

obj:
	mkdir -p obj

data:
	mkdir -p data

mlev-hiv-sim: $(OBJ_FILES)
	$(CC) $(LD_FLAGS) $^ $(LD_LIBS) -o $@

obj/%.o: src/%.cpp $^
	$(CC) $(CC_FLAGS) $< -o $@

CC_FLAGS += -MMD
-include $(OBJ_FILES:.o=.d)
