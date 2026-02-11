#
CUDA_BASE_DIR := /usr/local/cuda
GPU_ARCH := sm_89
ROOT_DIR := $(shell pwd)

# 默认不启用GPU（如果未指定）
DDPCA_USE_GPU ?= 0

# 编译器设置
LIB_DIR := $(ROOT_DIR)/lib

# 从基础路径派生include和lib路径 -ljemalloc --param=destructive-interference-size=64 -g -lpthread -lm -ldl -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-unused-variable $(ROOT_DIR):
ifeq ($(DDPCA_USE_GPU), 1)
	# GPU模式 - 使用nvcc编译器
	CXX = nvcc
	CXXFLAGS = -std=c++20 -O3 -arch=$(GPU_ARCH) -DNDEBUG \
               -Xcompiler -march=native -Xcompiler -pthread \
			   -Xcompiler -Wall -Xcompiler -Wextra \
               -DDDPCA_USE_GPU=$(DDPCA_USE_GPU) -I$(CUDA_BASE_DIR)/include
	LDFLAGS = -L$(LIB_DIR):$(CUDA_BASE_DIR)/lib64 -lDdpca -lcusparse
else
	# CPU模式 - 使用g++编译器
	CXX = g++
	CXXFLAGS = -std=c++20 -O3 -DNDEBUG -march=native -pthread -Wall -Wextra \
	           -DDDPCA_USE_GPU=$(DDPCA_USE_GPU)
	LDFLAGS = -L$(LIB_DIR) -lDdpca
endif
AR = ar
ARFLAGS = rcs

# 目录设置
CONTACT_DIR = Contact
DECOMPOSITION_DIR = Decomposition
GENERAL_DIR = General
MESH_DIR = Mesh
SOLVER_DIR = Solver
OBJ_DIR = obj

# 自动搜索所有.cpp文件
SRC_DIRS = $(CONTACT_DIR) $(DECOMPOSITION_DIR) $(GENERAL_DIR) $(MESH_DIR) $(SOLVER_DIR) $(SOLVER_DIR)/CSparse
SRCS := $(foreach dir,$(SRC_DIRS),$(shell find $(dir) -name "*.cpp"))
# 将源文件路径转换为对象文件路径
OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(notdir $(SRCS)))
# 静态库名称
LIBRARY = $(LIB_DIR)/libDdpca.a

export CXX CXXFLAGS LDFLAGS LIB_DIR

# 默认目标
all: $(LIBRARY) examples

# 创建静态库
$(LIBRARY): $(OBJS)
	@mkdir -p $(LIB_DIR)
	$(AR) $(ARFLAGS) $@ $^

# 编译规则
vpath %.cpp $(SRC_DIRS)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

examples: $(LIBRARY)
	@echo "Building examples..."
	$(MAKE) -C Examples

# 清理规则
clean: clean_examples
	@echo "Cleaning up root project..."
	rm -rf $(OBJ_DIR) $(LIB_DIR)

clean_examples:
	@$(MAKE) -C Examples clean

# 伪目标
.PHONY: all clean
