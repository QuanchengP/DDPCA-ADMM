
SRCS := Test.cpp
TARGET := $(patsubst %.cpp,%,$(SRCS))

EXAMPLE_DIR := examples

.PHONY: all clean

all:
	$(MAKE) -C $(EXAMPLE_DIR)
	g++ -O3 -std=c++17 -fopenmp -march=native $(SRCS) -o $(TARGET) -I./

clean:
	$(MAKE) -C $(EXAMPLE_DIR) clean
	rm -rf $(TARGET)
	@for dir in Beam Block Cylinder Dehw Torsion; do \
		if [ -d "$$dir" ]; then \
			echo "rm -rf '$$dir'..."; \
			rm -rf "$$dir"; \
		fi; \
	done
