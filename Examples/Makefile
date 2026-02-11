# 自动检测所有.cpp文件
CPP_SOURCES = $(wildcard *.cpp)
# 从.cpp文件生成对应的可执行文件名（去掉.cpp后缀）
EXECUTABLES = $(CPP_SOURCES:.cpp=)

# 默认目标：构建所有可执行文件
all: $(EXECUTABLES)

# 为每个可执行文件创建规则
%: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# 清理生成的可执行文件和同名+下划线的文件夹
clean:
	@echo "Cleaning up examples..."
	rm -f $(EXECUTABLES)
	rm -f testLogg.txt
	@for exe in $(EXECUTABLES); do \
		if [ -d "$${exe}_" ]; then \
			rm -rf "$${exe}_"; \
			echo "Removed directory: $${exe}_"; \
		fi \
	done

# 声明伪目标
.PHONY: all clean