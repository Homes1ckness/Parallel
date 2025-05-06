#include "md5.h"
#include <iomanip>
#include <assert.h>
#include <chrono>
#include <arm_neon.h>

using namespace std;
using namespace chrono;

/**
 * StringProcess: 将单个输入字符串转换成MD5计算所需的消息数组
 * @param input 输入
 * @param[out] n_byte 用于给调用者传递额外的返回值，即最终Byte数组的长度
 * @return Byte消息数组
 */
Byte *StringProcess(string input, int *n_byte)
{
	// 将输入的字符串转换为Byte为单位的数组
	Byte *blocks = (Byte *)input.c_str();
	int length = input.length();

	// 计算原始消息长度（以比特为单位）
	int bitLength = length * 8;

	// paddingBits: 原始消息需要的padding长度（以bit为单位）
	// 对于给定的消息，将其补齐至length%512==448为止
	// 需要注意的是，即便给定的消息满足length%512==448，也需要再pad 512bits
	int paddingBits = bitLength % 512;
	if (paddingBits > 448)
	{
		paddingBits = 512 - (paddingBits - 448);
	}
	else if (paddingBits < 448)
	{
		paddingBits = 448 - paddingBits;
	}
	else if (paddingBits == 448)
	{
		paddingBits = 512;
	}

	// 原始消息需要的padding长度（以Byte为单位）
	int paddingBytes = paddingBits / 8;
	// 创建最终的字节数组
	// length + paddingBytes + 8:
	// 1. length为原始消息的长度（bits）
	// 2. paddingBytes为原始消息需要的padding长度（Bytes）
	// 3. 在pad到length%512==448之后，需要额外附加64bits的原始消息长度，即8个bytes
	int paddedLength = length + paddingBytes + 8;
	Byte *paddedMessage = new Byte[paddedLength];

	// 复制原始消息
	memcpy(paddedMessage, blocks, length);

	// 添加填充字节。填充时，第一位为1，后面的所有位均为0。
	// 所以第一个byte是0x80
	paddedMessage[length] = 0x80;							 // 添加一个0x80字节
	memset(paddedMessage + length + 1, 0, paddingBytes - 1); // 填充0字节

	// 添加消息长度（64比特，小端格式）
	for (int i = 0; i < 8; ++i)
	{
		// 特别注意此处应当将bitLength转换为uint64_t
		// 这里的length是原始消息的长度
		paddedMessage[length + paddingBytes + i] = ((uint64_t)length * 8 >> (i * 8)) & 0xFF;
	}

	// 验证长度是否满足要求。此时长度应当是512bit的倍数
	int residual = 8 * paddedLength % 512;
	// assert(residual == 0);

	// 在填充+添加长度之后，消息被分为n_blocks个512bit的部分
	*n_byte = paddedLength;
	return paddedMessage;
}


/**
 * MD5Hash: 将单个输入字符串转换成MD5
 * @param input 输入
 * @param[out] state 用于给调用者传递额外的返回值，即最终的缓冲区，也就是MD5的结果
 * @return Byte消息数组
 */
void MD5Hash(string input, bit32 *state)
{

	Byte *paddedMessage;
	int *messageLength = new int[1];
	for (int i = 0; i < 1; i += 1)
	{
		paddedMessage = StringProcess(input, &messageLength[i]);
		// cout<<messageLength[i]<<endl;
		assert(messageLength[i] == messageLength[0]);
	}
	int n_blocks = messageLength[0] / 64;

	// bit32* state= new bit32[4];
	state[0] = 0x67452301;
	state[1] = 0xefcdab89;
	state[2] = 0x98badcfe;
	state[3] = 0x10325476;

	// 逐block地更新state
	for (int i = 0; i < n_blocks; i += 1)
	{
		bit32 x[16];

		// 下面的处理，在理解上较为复杂
		for (int i1 = 0; i1 < 16; ++i1)
		{
			x[i1] = (paddedMessage[4 * i1 + i * 64]) |
					(paddedMessage[4 * i1 + 1 + i * 64] << 8) |
					(paddedMessage[4 * i1 + 2 + i * 64] << 16) |
					(paddedMessage[4 * i1 + 3 + i * 64] << 24);
		}

		bit32 a = state[0], b = state[1], c = state[2], d = state[3];

		auto start = system_clock::now();
		/* Round 1 */
		FF(a, b, c, d, x[0], s11, 0xd76aa478);
		FF(d, a, b, c, x[1], s12, 0xe8c7b756);
		FF(c, d, a, b, x[2], s13, 0x242070db);
		FF(b, c, d, a, x[3], s14, 0xc1bdceee);
		FF(a, b, c, d, x[4], s11, 0xf57c0faf);
		FF(d, a, b, c, x[5], s12, 0x4787c62a);
		FF(c, d, a, b, x[6], s13, 0xa8304613);
		FF(b, c, d, a, x[7], s14, 0xfd469501);
		FF(a, b, c, d, x[8], s11, 0x698098d8);
		FF(d, a, b, c, x[9], s12, 0x8b44f7af);
		FF(c, d, a, b, x[10], s13, 0xffff5bb1);
		FF(b, c, d, a, x[11], s14, 0x895cd7be);
		FF(a, b, c, d, x[12], s11, 0x6b901122);
		FF(d, a, b, c, x[13], s12, 0xfd987193);
		FF(c, d, a, b, x[14], s13, 0xa679438e);
		FF(b, c, d, a, x[15], s14, 0x49b40821);

		/* Round 2 */
		GG(a, b, c, d, x[1], s21, 0xf61e2562);
		GG(d, a, b, c, x[6], s22, 0xc040b340);
		GG(c, d, a, b, x[11], s23, 0x265e5a51);
		GG(b, c, d, a, x[0], s24, 0xe9b6c7aa);
		GG(a, b, c, d, x[5], s21, 0xd62f105d);
		GG(d, a, b, c, x[10], s22, 0x2441453);
		GG(c, d, a, b, x[15], s23, 0xd8a1e681);
		GG(b, c, d, a, x[4], s24, 0xe7d3fbc8);
		GG(a, b, c, d, x[9], s21, 0x21e1cde6);
		GG(d, a, b, c, x[14], s22, 0xc33707d6);
		GG(c, d, a, b, x[3], s23, 0xf4d50d87);
		GG(b, c, d, a, x[8], s24, 0x455a14ed);
		GG(a, b, c, d, x[13], s21, 0xa9e3e905);
		GG(d, a, b, c, x[2], s22, 0xfcefa3f8);
		GG(c, d, a, b, x[7], s23, 0x676f02d9);
		GG(b, c, d, a, x[12], s24, 0x8d2a4c8a);

		/* Round 3 */
		HH(a, b, c, d, x[5], s31, 0xfffa3942);
		HH(d, a, b, c, x[8], s32, 0x8771f681);
		HH(c, d, a, b, x[11], s33, 0x6d9d6122);
		HH(b, c, d, a, x[14], s34, 0xfde5380c);
		HH(a, b, c, d, x[1], s31, 0xa4beea44);
		HH(d, a, b, c, x[4], s32, 0x4bdecfa9);
		HH(c, d, a, b, x[7], s33, 0xf6bb4b60);
		HH(b, c, d, a, x[10], s34, 0xbebfbc70);
		HH(a, b, c, d, x[13], s31, 0x289b7ec6);
		HH(d, a, b, c, x[0], s32, 0xeaa127fa);
		HH(c, d, a, b, x[3], s33, 0xd4ef3085);
		HH(b, c, d, a, x[6], s34, 0x4881d05);
		HH(a, b, c, d, x[9], s31, 0xd9d4d039);
		HH(d, a, b, c, x[12], s32, 0xe6db99e5);
		HH(c, d, a, b, x[15], s33, 0x1fa27cf8);
		HH(b, c, d, a, x[2], s34, 0xc4ac5665);

		/* Round 4 */
		II(a, b, c, d, x[0], s41, 0xf4292244);
		II(d, a, b, c, x[7], s42, 0x432aff97);
		II(c, d, a, b, x[14], s43, 0xab9423a7);
		II(b, c, d, a, x[5], s44, 0xfc93a039);
		II(a, b, c, d, x[12], s41, 0x655b59c3);
		II(d, a, b, c, x[3], s42, 0x8f0ccc92);
		II(c, d, a, b, x[10], s43, 0xffeff47d);
		II(b, c, d, a, x[1], s44, 0x85845dd1);
		II(a, b, c, d, x[8], s41, 0x6fa87e4f);
		II(d, a, b, c, x[15], s42, 0xfe2ce6e0);
		II(c, d, a, b, x[6], s43, 0xa3014314);
		II(b, c, d, a, x[13], s44, 0x4e0811a1);
		II(a, b, c, d, x[4], s41, 0xf7537e82);
		II(d, a, b, c, x[11], s42, 0xbd3af235);
		II(c, d, a, b, x[2], s43, 0x2ad7d2bb);
		II(b, c, d, a, x[9], s44, 0xeb86d391);

		state[0] += a;
		state[1] += b;
		state[2] += c;
		state[3] += d;
	}

	// 下面的处理，在理解上较为复杂
	for (int i = 0; i < 4; i++)
	{
		uint32_t value = state[i];
		state[i] = ((value & 0xff) << 24) |		 // 将最低字节移到最高位
				   ((value & 0xff00) << 8) |	 // 将次低字节左移
				   ((value & 0xff0000) >> 8) |	 // 将次高字节右移
				   ((value & 0xff000000) >> 24); // 将最高字节移到最低位
	}

	// 输出最终的hash结果
	// for (int i1 = 0; i1 < 4; i1 += 1)
	// {
	// 	cout << std::setw(8) << std::setfill('0') << hex << state[i1];
	// }
	// cout << endl;

	// 释放动态分配的内存
	// 实现SIMD并行算法的时候，也请记得及时回收内存！
	delete[] paddedMessage;
	delete[] messageLength;
}


// NEON版本的基本MD5函数
inline uint32x4_t F_SIMD(uint32x4_t x, uint32x4_t y, uint32x4_t z){
	uint32x4_t temp1 = vandq_u32(x,y); 		// x & y（按位与）
	uint32x4_t temp2 = vbicq_u32(z,x);    	// z & ~x（按位清除，即先对x取反，再与z按位与）
	return vorrq_u32(temp1,temp2);			// temp1 | temp2（按位或）
}

inline uint32x4_t G_SIMD(uint32x4_t x, uint32x4_t y, uint32x4_t z){
	uint32x4_t temp1 = vandq_u32(x,z);		// x & z
	uint32x4_t temp2 = vbicq_u32(y,z);		// y & ~z
	return vorrq_u32(temp1,temp2);			// temp1 | temp2
}

inline uint32x4_t H_SIMD(uint32x4_t x, uint32x4_t y, uint32x4_t z){
	return veorq_u32(veorq_u32(x, y), z);		//x ^ y ^ z  异或
}

inline uint32x4_t I_SIMD(uint32x4_t x, uint32x4_t y, uint32x4_t z){
	// uint32x4_t temp = vbicq_u32(x, z);
	uint32x4_t temp = vorrq_u32(x, vmvnq_u32(z));  // x | (~z)
	return veorq_u32(y, temp);
}

//NEON版本的左循环移位函数
// 为什么vshlq_n_u32和vshrq_n_u32会报错??
inline uint32x4_t ROTATELEFT_SIMD(uint32x4_t num, int n) {

    return vorrq_u32(vshlq_n_u32(num, n), vshrq_n_u32(num, (32-n)));
}


// NEON版本的四轮运算函数
// vdupq_n_u32(ac):将一个 32 位无符号整型标量 ac 复制到一个 uint32x4_t 类型向量 的所有 4 个元素中。
inline void FF_SIMD(uint32x4_t &a, uint32x4_t b, uint32x4_t c, uint32x4_t d, uint32x4_t x, int s, uint32_t ac){
	uint32x4_t ac_vec = vdupq_n_u32(ac);
	a = vaddq_u32(a, vaddq_u32(vaddq_u32(F_SIMD(b, c, d), x), ac_vec)); //(a) += F ((b), (c), (d)) + (x) + ac;
	a = ROTATELEFT_SIMD(a, s);
	a = vaddq_u32(a, b);
}

inline void GG_SIMD(uint32x4_t &a, uint32x4_t b, uint32x4_t c, uint32x4_t d, uint32x4_t x, int s, uint32_t ac) {
    uint32x4_t ac_vec = vdupq_n_u32(ac);
    a = vaddq_u32(a, vaddq_u32(vaddq_u32(G_SIMD(b, c, d), x), ac_vec));
    a = ROTATELEFT_SIMD(a, s);
    a = vaddq_u32(a, b);
}

inline void HH_SIMD(uint32x4_t &a, uint32x4_t b, uint32x4_t c, uint32x4_t d, uint32x4_t x, int s, uint32_t ac) {
    uint32x4_t ac_vec = vdupq_n_u32(ac);
    a = vaddq_u32(a, vaddq_u32(vaddq_u32(H_SIMD(b, c, d), x), ac_vec));
    a = ROTATELEFT_SIMD(a, s);
    a = vaddq_u32(a, b);
}

inline void II_SIMD(uint32x4_t &a, uint32x4_t b, uint32x4_t c, uint32x4_t d, uint32x4_t x, int s, uint32_t ac) {
    uint32x4_t ac_vec = vdupq_n_u32(ac);
    a = vaddq_u32(a, vaddq_u32(vaddq_u32(I_SIMD(b, c, d), x), ac_vec));
    a = ROTATELEFT_SIMD(a, s);
    a = vaddq_u32(a, b);
}


// SIMD版本的MD5哈希函数
void MD5Hash_SIMD(vector<string> inputs, vector<bit32*> states){
	// 同时处理4个消息（NEON的uint32x4_t可以同时处理4个32位整数）
	const int batchSize = 4;

	if (inputs.size() < batchSize){
		// 如果输入少于4个，就按照普通哈希处理
		for (size_t i = 0; i < inputs.size();i++){
			MD5Hash(inputs[i], states[i]);
		}
		return;
	}

	// 处理每一批4个消息
	for (size_t batchStart = 0; batchStart < inputs.size(); batchStart += batchSize) {
        // 确保不会越界
        if (batchStart + batchSize > inputs.size()) {
            // 处理剩余的消息
            for (size_t i = batchStart; i < inputs.size(); i++) {
                MD5Hash(inputs[i], states[i]);
            }
            break;
        }

		// 预处理消息
        vector<Byte*> paddedMessages;
        vector<int> messageLengths;
        
        for (size_t i = 0; i < batchSize; i++) {
            int length;
            Byte* paddedMessage = StringProcess(inputs[batchStart + i], &length);
            paddedMessages.push_back(paddedMessage);
            messageLengths.push_back(length);
        }

		 // 检查所有消息长度是否相同
		 bool sameLengths = true;
		 for (size_t i = 1; i < batchSize; i++) {
			 if (messageLengths[i] != messageLengths[0]) {
				 sameLengths = false;
				 break;
			 }
		 }

		 if (!sameLengths) {
            // 如果长度不同，使用普通方法处理
            for (size_t i = 0; i < batchSize; i++) {
                MD5Hash(inputs[batchStart + i], states[batchStart + i]);
                delete[] paddedMessages[i];
            }
            continue;
        }

		// 初始化MD5状态
        for (size_t i = 0; i < batchSize; i++) {
            states[batchStart + i][0] = 0x67452301;
            states[batchStart + i][1] = 0xefcdab89;
            states[batchStart + i][2] = 0x98badcfe;
            states[batchStart + i][3] = 0x10325476;
        }

		int n_blocks = messageLengths[0] / 64;

		// 对每个块进行处理
        for (int blockIdx = 0; blockIdx < n_blocks; blockIdx++) {
            // 准备SIMD数据
            uint32x4_t x[16];
            
			for (int i = 0; i < 16; i++) {
				alignas(16) uint32_t values[4];  // 16字节对齐
				for (int j = 0; j < 4; j++) {
					int msgIdx = j;
					int offset = 4 * i + blockIdx * 64;
					values[j] = (paddedMessages[j][offset]) |
							   (paddedMessages[j][offset + 1] << 8) |
							   (paddedMessages[j][offset + 2] << 16) |
							   (paddedMessages[j][offset + 3] << 24);
				}
				x[i] = vld1q_u32(values);  // 直接使用NEON加载指令
			}
			// 加载当前状态
			alignas(16) uint32_t a_values[4], b_values[4], c_values[4], d_values[4];
			for (int j = 0; j < 4; j++) {
				a_values[j] = states[batchStart + j][0];
				b_values[j] = states[batchStart + j][1];
				c_values[j] = states[batchStart + j][2];
				d_values[j] = states[batchStart + j][3];
			}
			
            
            uint32x4_t a = vld1q_u32(a_values);
            uint32x4_t b = vld1q_u32(b_values);
            uint32x4_t c = vld1q_u32(c_values);
            uint32x4_t d = vld1q_u32(d_values);

			// uint32x4_t a, b, c, d;
			// memcpy(&a, a_values, sizeof(a_values));
			// memcpy(&b, b_values, sizeof(b_values));
			// memcpy(&c, c_values, sizeof(c_values));
			// memcpy(&d, d_values, sizeof(d_values));
            
            // 保存原始状态
            uint32x4_t original_a = a;
            uint32x4_t original_b = b;
            uint32x4_t original_c = c;
            uint32x4_t original_d = d;

			 // MD5四轮变换
            /* Round 1 */
            FF_SIMD(a, b, c, d, x[0], s11, 0xd76aa478);
            FF_SIMD(d, a, b, c, x[1], s12, 0xe8c7b756);
            FF_SIMD(c, d, a, b, x[2], s13, 0x242070db);
            FF_SIMD(b, c, d, a, x[3], s14, 0xc1bdceee);
            FF_SIMD(a, b, c, d, x[4], s11, 0xf57c0faf);
            FF_SIMD(d, a, b, c, x[5], s12, 0x4787c62a);
            FF_SIMD(c, d, a, b, x[6], s13, 0xa8304613);
            FF_SIMD(b, c, d, a, x[7], s14, 0xfd469501);
            FF_SIMD(a, b, c, d, x[8], s11, 0x698098d8);
            FF_SIMD(d, a, b, c, x[9], s12, 0x8b44f7af);
            FF_SIMD(c, d, a, b, x[10], s13, 0xffff5bb1);
            FF_SIMD(b, c, d, a, x[11], s14, 0x895cd7be);
            FF_SIMD(a, b, c, d, x[12], s11, 0x6b901122);
            FF_SIMD(d, a, b, c, x[13], s12, 0xfd987193);
            FF_SIMD(c, d, a, b, x[14], s13, 0xa679438e);
            FF_SIMD(b, c, d, a, x[15], s14, 0x49b40821);

            /* Round 2 */
            GG_SIMD(a, b, c, d, x[1], s21, 0xf61e2562);
            GG_SIMD(d, a, b, c, x[6], s22, 0xc040b340);
            GG_SIMD(c, d, a, b, x[11], s23, 0x265e5a51);
            GG_SIMD(b, c, d, a, x[0], s24, 0xe9b6c7aa);
            GG_SIMD(a, b, c, d, x[5], s21, 0xd62f105d);
            GG_SIMD(d, a, b, c, x[10], s22, 0x2441453);
            GG_SIMD(c, d, a, b, x[15], s23, 0xd8a1e681);
            GG_SIMD(b, c, d, a, x[4], s24, 0xe7d3fbc8);
            GG_SIMD(a, b, c, d, x[9], s21, 0x21e1cde6);
            GG_SIMD(d, a, b, c, x[14], s22, 0xc33707d6);
            GG_SIMD(c, d, a, b, x[3], s23, 0xf4d50d87);
            GG_SIMD(b, c, d, a, x[8], s24, 0x455a14ed);
            GG_SIMD(a, b, c, d, x[13], s21, 0xa9e3e905);
            GG_SIMD(d, a, b, c, x[2], s22, 0xfcefa3f8);
            GG_SIMD(c, d, a, b, x[7], s23, 0x676f02d9);
            GG_SIMD(b, c, d, a, x[12], s24, 0x8d2a4c8a);

            /* Round 3 */
            HH_SIMD(a, b, c, d, x[5], s31, 0xfffa3942);
            HH_SIMD(d, a, b, c, x[8], s32, 0x8771f681);
            HH_SIMD(c, d, a, b, x[11], s33, 0x6d9d6122);
            HH_SIMD(b, c, d, a, x[14], s34, 0xfde5380c);
            HH_SIMD(a, b, c, d, x[1], s31, 0xa4beea44);
            HH_SIMD(d, a, b, c, x[4], s32, 0x4bdecfa9);
            HH_SIMD(c, d, a, b, x[7], s33, 0xf6bb4b60);
            HH_SIMD(b, c, d, a, x[10], s34, 0xbebfbc70);
            HH_SIMD(a, b, c, d, x[13], s31, 0x289b7ec6);
            HH_SIMD(d, a, b, c, x[0], s32, 0xeaa127fa);
            HH_SIMD(c, d, a, b, x[3], s33, 0xd4ef3085);
            HH_SIMD(b, c, d, a, x[6], s34, 0x4881d05);
            HH_SIMD(a, b, c, d, x[9], s31, 0xd9d4d039);
            HH_SIMD(d, a, b, c, x[12], s32, 0xe6db99e5);
            HH_SIMD(c, d, a, b, x[15], s33, 0x1fa27cf8);
            HH_SIMD(b, c, d, a, x[2], s34, 0xc4ac5665);

            /* Round 4 */
            II_SIMD(a, b, c, d, x[0], s41, 0xf4292244);
            II_SIMD(d, a, b, c, x[7], s42, 0x432aff97);
            II_SIMD(c, d, a, b, x[14], s43, 0xab9423a7);
            II_SIMD(b, c, d, a, x[5], s44, 0xfc93a039);
            II_SIMD(a, b, c, d, x[12], s41, 0x655b59c3);
            II_SIMD(d, a, b, c, x[3], s42, 0x8f0ccc92);
            II_SIMD(c, d, a, b, x[10], s43, 0xffeff47d);
            II_SIMD(b, c, d, a, x[1], s44, 0x85845dd1);
            II_SIMD(a, b, c, d, x[8], s41, 0x6fa87e4f);
            II_SIMD(d, a, b, c, x[15], s42, 0xfe2ce6e0);
            II_SIMD(c, d, a, b, x[6], s43, 0xa3014314);
            II_SIMD(b, c, d, a, x[13], s44, 0x4e0811a1);
            II_SIMD(a, b, c, d, x[4], s41, 0xf7537e82);
            II_SIMD(d, a, b, c, x[11], s42, 0xbd3af235);
            II_SIMD(c, d, a, b, x[2], s43, 0x2ad7d2bb);
            II_SIMD(b, c, d, a, x[9], s44, 0xeb86d391);
            
            // 更新状态
            a = vaddq_u32(a, original_a);
            b = vaddq_u32(b, original_b);
            c = vaddq_u32(c, original_c);
            d = vaddq_u32(d, original_d);

			// 存储结果
            vst1q_u32(a_values, a);
            vst1q_u32(b_values, b);
            vst1q_u32(c_values, c);
            vst1q_u32(d_values, d);

			for (int j = 0; j < 4; j++) {
                states[batchStart + j][0] = a_values[j];
                states[batchStart + j][1] = b_values[j];
                states[batchStart + j][2] = c_values[j];
                states[batchStart + j][3] = d_values[j];
            }
        }

		// 字节序调整
        for (size_t i = 0; i < batchSize; i++) {
            for (int j = 0; j < 4; j++) {
                uint32_t value = states[batchStart + i][j];
                states[batchStart + i][j] = ((value & 0xff) << 24) |
                                          ((value & 0xff00) << 8) |
                                          ((value & 0xff0000) >> 8) |
                                          ((value & 0xff000000) >> 24);
            }
        }

		// 释放内存
        for (size_t i = 0; i < batchSize; i++) {
            delete[] paddedMessages[i];
        }
    }

}
