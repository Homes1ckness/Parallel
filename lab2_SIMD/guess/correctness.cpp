#include "PCFG.h"
#include <chrono>
#include <fstream>
#include "md5.h"
#include <iomanip>
using namespace std;
using namespace chrono;

// 编译指令如下：
// g++ correctness.cpp train.cpp guessing.cpp md5.cpp -o test

// 通过这个函数，你可以验证你实现的SIMD哈希函数的正确性
int main()
{
    string testStr = "bvaisdbjasdkafkasdfnavkjnakdjfejfanjsdnfkajdfkajdfjkwanfdjaknsvjkanbjbjadfajwefajksdfakdnsvjadfasjdvabvaisdbjasdkafkasdfnavkjnakdjfejfanjsdnfkajdfkajdfjkwanfdjaknsvjkanbjbjadfajwefajksdfakdnsvjadfasjdvabvaisdbjasdkafkasdfnavkjnakdjfejfanjsdnfkajdfkajdfjkwanfdjaknsvjkanbjbjadfajwefajksdfakdnsvjadfasjdvabvaisdbjasdkafkasdfnavkjnakdjfejfanjsdnfkajdfkajdfjkwanfdjaknsvjkanbjbjadfajwefajksdfakdnsvjadfasjdva";
    
    // 标准MD5计算
    bit32 state[4];
    MD5Hash(testStr, state);
    
    cout << "标准MD5结果: ";
    for (int i1 = 0; i1 < 4; i1 += 1)
    {
        cout << std::setw(8) << std::setfill('0') << hex << state[i1];
    }
    cout << endl;
    
    // SIMD版本MD5计算
    vector<string> inputs(4, testStr); // 4个相同的输入
    vector<bit32*> simdStates;
    for (int i = 0; i < 4; i++) {
        simdStates.push_back(new bit32[4]);
    }
    
    MD5Hash_SIMD(inputs, simdStates);
    
    cout << "SIMD MD5结果: ";
    for (int i1 = 0; i1 < 4; i1 += 1)
    {
        cout << std::setw(8) << std::setfill('0') << hex << simdStates[0][i1];
    }
    cout << endl;
    
    // 检查结果是否一致
    bool match = true;
    for (int i = 0; i < 4; i++) {
        if (state[i] != simdStates[0][i]) {
            match = false;
            break;
        }
    }
    
    cout << "结果比较: " << (match ? "一致 ✓" : "不一致 ✗") << endl;
}