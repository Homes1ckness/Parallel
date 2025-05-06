#include "PCFG.h"
#include <chrono>
#include <fstream>
#include "md5.h"
#include <iomanip>
using namespace std;
using namespace chrono;

// 编译指令如下
// g++ main.cpp train.cpp guessing.cpp md5.cpp -o main
// g++ main.cpp train.cpp guessing.cpp md5.cpp -o main -O1
// g++ main.cpp train.cpp guessing.cpp md5.cpp -o main -O2

int main()
{
    double time_hash = 0;  // 用于MD5哈希的时间
    double time_guess = 0; // 哈希和猜测的总时长
    double time_train = 0; // 模型训练的总时长
    PriorityQueue q;
    auto start_train = system_clock::now();
    q.m.train("/guessdata/Rockyou-singleLined-full.txt");
    q.m.order();
    auto end_train = system_clock::now();
    auto duration_train = duration_cast<microseconds>(end_train - start_train);
    time_train = double(duration_train.count()) * microseconds::period::num / microseconds::period::den;

    q.init();
    cout << "here" << endl;
    int curr_num = 0;
    auto start = system_clock::now();
    // 由于需要定期清空内存，我们在这里记录已生成的猜测总数
    int history = 0;
    // std::ofstream a("./files/results.txt");
    while (!q.priority.empty())
    {
        q.PopNext();
        q.total_guesses = q.guesses.size();
        if (q.total_guesses - curr_num >= 100000)
        {
            cout << "Guesses generated: " <<history + q.total_guesses << endl;
            curr_num = q.total_guesses;

            // 在此处更改实验生成的猜测上限
            int generate_n=10000000;
            if (history + q.total_guesses > 10000000)
            {
                auto end = system_clock::now();
                auto duration = duration_cast<microseconds>(end - start);
                time_guess = double(duration.count()) * microseconds::period::num / microseconds::period::den;
                cout << "Guess time:" << time_guess - time_hash << "seconds"<< endl;
                cout << "Hash time:" << time_hash << "seconds"<<endl;
                cout << "Train time:" << time_train <<"seconds"<<endl;
                break;
            }
        }
        // 为了避免内存超限，我们在q.guesses中口令达到一定数目时，将其中的所有口令取出并且进行哈希
        // 然后，q.guesses将会被清空。为了有效记录已经生成的口令总数，维护一个history变量来进行记录
        if (curr_num > 1000000)
        {
            auto start_hash = system_clock::now();
            // bit32 state[4];
            // for (string pw : q.guesses)
            // {
            //     // TODO：对于SIMD实验，将这里替换成你的SIMD MD5函数
            //     MD5Hash(pw, state);

            //     // 以下注释部分用于输出猜测和哈希，但是由于自动测试系统不太能写文件，所以这里你可以改成cout
            //     // a<<pw<<"\t";
            //     // for (int i1 = 0; i1 < 4; i1 += 1)
            //     // {
            //     //     a << std::setw(8) << std::setfill('0') << hex << state[i1];
            //     // }
            //     // a << endl;
            // }

            // SIMD版本 - 批量处理密码
            const int batchSize = 4; // NEON处理4个32位整数
            int numPasswords = q.guesses.size();

            
            
            for (int i = 0; i < numPasswords; i += batchSize) {
                // 确定当前批次大小(避免越界)
                int currentBatchSize = min(batchSize, numPasswords - i);
                
                // 准备当前批次的输入和状态数组
                vector<string> batchInputs;
                vector<bit32*> batchStates(currentBatchSize);
                
                for (int j = 0; j < currentBatchSize; j++) {
                    batchInputs.push_back(q.guesses[i + j]);
                    batchStates[j] = new bit32[4];
                }

               // 修改条件：只要有至少1个密码就可以处理
                if (currentBatchSize >= 1) {
                    // 如果不是4的倍数，添加空字符串填充到4个
                    while (batchInputs.size() < batchSize) {
                        batchInputs.push_back("");
                        batchStates.push_back(new bit32[4]);
                    }
                    MD5Hash_SIMD(batchInputs, batchStates);
                }
                
                // 释放内存
                for (size_t j = 0; j < batchStates.size(); j++) {
                    delete[] batchStates[j];
                }
            }

            // 在这里对哈希所需的总时长进行计算
            auto end_hash = system_clock::now();
            auto duration = duration_cast<microseconds>(end_hash - start_hash);
            time_hash += double(duration.count()) * microseconds::period::num / microseconds::period::den;

            // 记录已经生成的口令总数
            history += curr_num;
            curr_num = 0;
            q.guesses.clear();
        }
    }
}
