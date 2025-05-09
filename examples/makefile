
SRCS := $(wildcard *.cpp)

TARGET := $(SRCS:.cpp=)

all: $(TARGET)

%: %.cpp
	g++ -O3 -std=c++17 -fopenmp -march=native $< -o $@ -I../

clean:
	rm -f $(TARGET)

.PHONY: all clean
