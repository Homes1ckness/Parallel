# GPU并行化的PCFG口令猜测程序

## 实现说明

本程序对原有的两个关键for循环进行了GPU并行化改造：

### 1. 单个segment的情况（第一个循环）
- 原始代码：串行遍历所有values，逐个加入guesses
- GPU并行化：使用CUDA kernel `generate_single_segment_kernel` 在GPU上并行处理所有values

### 2. 多个segment的情况（第二个循环）  
- 原始代码：串行将prefix与每个value拼接，逐个加入guesses
- GPU并行化：使用CUDA kernel `generate_multi_segment_kernel` 在GPU上并行处理prefix和values的拼接

## 关键组件

### CUDA Kernels
1. `generate_single_segment_kernel`: 处理单个segment的并行复制
2. `generate_multi_segment_kernel`: 处理多个segment的并行字符串拼接

### 辅助函数
1. `copy_strings_to_device`: 将vector<string>复制到GPU内存
2. `gpu_generate_single_segment`: 单个segment的GPU生成函数
3. `gpu_generate_multi_segment`: 多个segment的GPU生成函数

## 编译方法

### 使用nvcc编译（推荐）
```bash
nvcc -x cu guessing.cpp train.cpp md5.cpp -o main_gpu -I. -lcudart
```

### 使用nvcc + g++混合编译
```bash
nvcc -c guessing.cpp -o guessing.o
g++ main.cpp train.cpp guessing.o md5.cpp -lcudart -o main_gpu
```

### 带优化的编译
```bash
nvcc -x cu guessing.cpp train.cpp md5.cpp -o main_gpu -I. -lcudart -O2
```

## 使用要求

1. **硬件要求**：需要CUDA兼容的NVIDIA GPU
2. **软件要求**：安装CUDA工具包(CUDA Toolkit)
3. **驱动要求**：安装匹配的NVIDIA驱动程序

## 性能优势

- **并行度**：将原本串行的字符串操作转换为高度并行的GPU计算
- **内存效率**：通过批量处理减少CPU-GPU数据传输开销
- **计算加速**：利用GPU的大量并行计算单元显著提升字符串生成速度

## 运行方法

```bash
./main_gpu
```

程序会自动使用GPU加速进行口令猜测生成。

## 注意事项

1. 如果没有GPU或CUDA环境，程序会报错
2. 确保GPU有足够的显存处理大量字符串数据
3. 程序会自动管理GPU内存的分配和释放
