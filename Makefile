# 编译器设置
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
TARGET = GeneralNewtonMethod
SOURCE = GeneralNewtonMethod.cpp

# 默认目标
all: $(TARGET)

# 编译目标
$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)

# 运行程序
run: $(TARGET)
	./$(TARGET)

# 调试版本
debug: CXXFLAGS += -g -DDEBUG
debug: $(TARGET)

# 清理生成的文件
clean:
	rm -f $(TARGET)

# 显示帮助信息
help:
	@echo "可用命令:"
	@echo "  make       - 编译程序"
	@echo "  make run   - 编译并运行程序"
	@echo "  make debug - 编译调试版本"
	@echo "  make clean - 清理生成的文件"
	@echo "  make help  - 显示此帮助信息"

# 防止与同名文件冲突
.PHONY: all run debug clean help