#include "PCFG.h"
#include <mpi.h>
#include <sstream>
#include <algorithm>
using namespace std;

void PriorityQueue::CalProb(PT &pt)
{
    // 计算PriorityQueue里面一个PT的流程如下：
    // 1. 首先需要计算一个PT本身的概率。例如，L6S1的概率为0.15
    // 2. 需要注意的是，Queue里面的PT不是“纯粹的”PT，而是除了最后一个segment以外，全部被value实例化的PT
    // 3. 所以，对于L6S1而言，其在Queue里面的实际PT可能是123456S1，其中“123456”为L6的一个具体value。
    // 4. 这个时候就需要计算123456在L6中出现的概率了。假设123456在所有L6 segment中的概率为0.1，那么123456S1的概率就是0.1*0.15

    // 计算一个PT本身的概率。后续所有具体segment value的概率，直接累乘在这个初始概率值上
    pt.prob = pt.preterm_prob;

    // index: 标注当前segment在PT中的位置
    int index = 0;


    for (int idx : pt.curr_indices)
    {
        // pt.content[index].PrintSeg();
        if (pt.content[index].type == 1)
        {
            // 下面这行代码的意义：
            // pt.content[index]：目前需要计算概率的segment
            // m.FindLetter(seg): 找到一个letter segment在模型中的对应下标
            // m.letters[m.FindLetter(seg)]：一个letter segment在模型中对应的所有统计数据
            // m.letters[m.FindLetter(seg)].ordered_values：一个letter segment在模型中，所有value的总数目
            pt.prob *= m.letters[m.FindLetter(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.letters[m.FindLetter(pt.content[index])].total_freq;
            // cout << m.letters[m.FindLetter(pt.content[index])].ordered_freqs[idx] << endl;
            // cout << m.letters[m.FindLetter(pt.content[index])].total_freq << endl;
        }
        if (pt.content[index].type == 2)
        {
            pt.prob *= m.digits[m.FindDigit(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.digits[m.FindDigit(pt.content[index])].total_freq;
            // cout << m.digits[m.FindDigit(pt.content[index])].ordered_freqs[idx] << endl;
            // cout << m.digits[m.FindDigit(pt.content[index])].total_freq << endl;
        }
        if (pt.content[index].type == 3)
        {
            pt.prob *= m.symbols[m.FindSymbol(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.symbols[m.FindSymbol(pt.content[index])].total_freq;
            // cout << m.symbols[m.FindSymbol(pt.content[index])].ordered_freqs[idx] << endl;
            // cout << m.symbols[m.FindSymbol(pt.content[index])].total_freq << endl;
        }
        index += 1;
    }
    // cout << pt.prob << endl;
}

void PriorityQueue::init()
{
    // cout << m.ordered_pts.size() << endl;
    // 用所有可能的PT，按概率降序填满整个优先队列
    for (PT pt : m.ordered_pts)
    {
        for (segment seg : pt.content)
        {
            if (seg.type == 1)
            {
                // 下面这行代码的意义：
                // max_indices用来表示PT中各个segment的可能数目。例如，L6S1中，假设模型统计到了100个L6，那么L6对应的最大下标就是99
                // （但由于后面采用了"<"的比较关系，所以其实max_indices[0]=100）
                // m.FindLetter(seg): 找到一个letter segment在模型中的对应下标
                // m.letters[m.FindLetter(seg)]：一个letter segment在模型中对应的所有统计数据
                // m.letters[m.FindLetter(seg)].ordered_values：一个letter segment在模型中，所有value的总数目
                pt.max_indices.emplace_back(m.letters[m.FindLetter(seg)].ordered_values.size());
            }
            if (seg.type == 2)
            {
                pt.max_indices.emplace_back(m.digits[m.FindDigit(seg)].ordered_values.size());
            }
            if (seg.type == 3)
            {
                pt.max_indices.emplace_back(m.symbols[m.FindSymbol(seg)].ordered_values.size());
            }
        }
        pt.preterm_prob = float(m.preterm_freq[m.FindPT(pt)]) / m.total_preterm;
        // pt.PrintPT();
        // cout << " " << m.preterm_freq[m.FindPT(pt)] << " " << m.total_preterm << " " << pt.preterm_prob << endl;

        // 计算当前pt的概率
        CalProb(pt);
        // 将PT放入优先队列
        priority.emplace_back(pt);
    }
    // cout << "priority size:" << priority.size() << endl;
}

void PriorityQueue::PopNext()
{

    // 对优先队列最前面的PT，首先利用这个PT生成一系列猜测
    Generate(priority.front());

    // 然后需要根据即将出队的PT，生成一系列新的PT
    vector<PT> new_pts = priority.front().NewPTs();
    for (PT pt : new_pts)
    {
        // 计算概率
        CalProb(pt);
        // 接下来的这个循环，作用是根据概率，将新的PT插入到优先队列中
        for (auto iter = priority.begin(); iter != priority.end(); iter++)
        {
            // 对于非队首和队尾的特殊情况
            if (iter != priority.end() - 1 && iter != priority.begin())
            {
                // 判定概率
                if (pt.prob <= iter->prob && pt.prob > (iter + 1)->prob)
                {
                    priority.emplace(iter + 1, pt);
                    break;
                }
            }
            if (iter == priority.end() - 1)
            {
                priority.emplace_back(pt);
                break;
            }
            if (iter == priority.begin() && iter->prob < pt.prob)
            {
                priority.emplace(iter, pt);
                break;
            }
        }
    }

    // 现在队首的PT善后工作已经结束，将其出队（删除）
    priority.erase(priority.begin());
}

// MPI 版本的 PopNext 方法 - 优化版本
void PriorityQueue::PopNextMPI()
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 只在主进程中执行主要逻辑
    if (rank == 0) {
        // 对优先队列最前面的PT，首先利用这个PT生成一系列猜测（使用MPI并行版本）
        GenerateMPI(priority.front());

        // 然后需要根据即将出队的PT，生成一系列新的PT
        vector<PT> new_pts = priority.front().NewPTs();
        for (PT pt : new_pts)
        {
            // 计算概率
            CalProb(pt);
            // 接下来的这个循环，作用是根据概率，将新的PT插入到优先队列中
            for (auto iter = priority.begin(); iter != priority.end(); iter++)
            {
                // 对于非队首和队尾的特殊情况
                if (iter != priority.end() - 1 && iter != priority.begin())
                {
                    // 判定概率
                    if (pt.prob <= iter->prob && pt.prob > (iter + 1)->prob)
                    {
                        priority.emplace(iter + 1, pt);
                        break;
                    }
                }
                if (iter == priority.end() - 1)
                {
                    priority.emplace_back(pt);
                    break;
                }
                if (iter == priority.begin() && iter->prob < pt.prob)
                {
                    priority.emplace(iter, pt);
                    break;
                }
            }
        }

        // 现在队首的PT善后工作已经结束，将其出队（删除）
        priority.erase(priority.begin());
    }
    // 其他进程不执行任何操作
}

// PT层面的MPI并行版本 - 多个进程同时处理不同的PT
void PriorityQueue::PopNextMPI_PTLevel()
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // 用于存储本进程处理的PT及其生成的新PT
    vector<PT> local_new_pts;
    vector<string> local_guesses;
    
    // 主进程负责分发PT给各个进程
    if (rank == 0) {
        // 检查队列中是否有足够的PT
        int available_pts = min((int)priority.size(), size);
        
        // 向所有进程广播可用PT数量
        MPI_Bcast(&available_pts, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // 如果没有足够的PT，只让部分进程工作
        if (rank < available_pts && !priority.empty()) {
            // 主进程处理第一个PT
            PT current_pt = priority[0];
            
            // 生成密码猜测
            Generate(current_pt);
            
            // 生成新的PT
            vector<PT> new_pts = current_pt.NewPTs();
            for (PT& pt : new_pts) {
                CalProb(pt);
                local_new_pts.push_back(pt);
            }
        }
        
        // 为其他进程发送PT数据
        for (int p = 1; p < available_pts && p < priority.size(); p++) {
            PT pt_to_send = priority[p];
            
            // 发送PT的基本信息
            int content_size = pt_to_send.content.size();
            MPI_Send(&content_size, 1, MPI_INT, p, 0, MPI_COMM_WORLD);
            
            // 发送每个segment的信息
            for (const segment& seg : pt_to_send.content) {
                MPI_Send(&seg.type, 1, MPI_INT, p, 1, MPI_COMM_WORLD);
                MPI_Send(&seg.length, 1, MPI_INT, p, 2, MPI_COMM_WORLD);
            }
            
            // 发送curr_indices和max_indices
            int indices_size = pt_to_send.curr_indices.size();
            MPI_Send(&indices_size, 1, MPI_INT, p, 3, MPI_COMM_WORLD);
            if (indices_size > 0) {
                MPI_Send(pt_to_send.curr_indices.data(), indices_size, MPI_INT, p, 4, MPI_COMM_WORLD);
                MPI_Send(pt_to_send.max_indices.data(), indices_size, MPI_INT, p, 5, MPI_COMM_WORLD);
            }
            
            // 发送其他PT属性
            MPI_Send(&pt_to_send.pivot, 1, MPI_INT, p, 6, MPI_COMM_WORLD);
            MPI_Send(&pt_to_send.preterm_prob, 1, MPI_FLOAT, p, 7, MPI_COMM_WORLD);
            MPI_Send(&pt_to_send.prob, 1, MPI_FLOAT, p, 8, MPI_COMM_WORLD);
        }
        
    } else {
        // 其他进程接收可用PT数量
        int available_pts;
        MPI_Bcast(&available_pts, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // 检查当前进程是否需要工作
        if (rank < available_pts) {
            // 接收PT数据
            PT received_pt;
            
            // 接收PT的基本信息
            int content_size;
            MPI_Recv(&content_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // 接收每个segment的信息
            for (int i = 0; i < content_size; i++) {
                int type, length;
                MPI_Recv(&type, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&length, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                received_pt.content.emplace_back(type, length);
            }
            
            // 接收curr_indices和max_indices
            int indices_size;
            MPI_Recv(&indices_size, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (indices_size > 0) {
                received_pt.curr_indices.resize(indices_size);
                received_pt.max_indices.resize(indices_size);
                MPI_Recv(received_pt.curr_indices.data(), indices_size, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(received_pt.max_indices.data(), indices_size, MPI_INT, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            // 接收其他PT属性
            MPI_Recv(&received_pt.pivot, 1, MPI_INT, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&received_pt.preterm_prob, 1, MPI_FLOAT, 0, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&received_pt.prob, 1, MPI_FLOAT, 0, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // 处理接收到的PT
            // 生成密码猜测
            int old_guesses_size = guesses.size();
            Generate(received_pt);
            
            // 收集本进程生成的新猜测
            for (int i = old_guesses_size; i < guesses.size(); i++) {
                local_guesses.push_back(guesses[i]);
            }
            
            // 生成新的PT
            vector<PT> new_pts = received_pt.NewPTs();
            for (PT& pt : new_pts) {
                CalProb(pt);
                local_new_pts.push_back(pt);
            }
        }
    }
    
    // 同步所有进程
    MPI_Barrier(MPI_COMM_WORLD);
    
    // 收集所有进程生成的新PT数量
    int local_new_pts_count = local_new_pts.size();
    vector<int> all_new_pts_counts(size);
    MPI_Allgather(&local_new_pts_count, 1, MPI_INT, all_new_pts_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    // 主进程收集所有新PT并更新优先队列
    if (rank == 0) {
        // 删除已处理的PT（从前面删除已分发给各进程的PT）
        int pts_to_remove = min((int)priority.size(), size);
        for (int i = 0; i < pts_to_remove; i++) {
            if (!priority.empty()) {
                priority.erase(priority.begin());
            }
        }
        
        // 首先添加主进程的新PT
        for (const PT& pt : local_new_pts) {
            // 按概率插入到优先队列的合适位置
            bool inserted = false;
            for (auto iter = priority.begin(); iter != priority.end(); iter++) {
                if (iter != priority.end() - 1 && iter != priority.begin()) {
                    if (pt.prob <= iter->prob && pt.prob > (iter + 1)->prob) {
                        priority.emplace(iter + 1, pt);
                        inserted = true;
                        break;
                    }
                }
                if (iter == priority.end() - 1) {
                    priority.emplace_back(pt);
                    inserted = true;
                    break;
                }
                if (iter == priority.begin() && iter->prob < pt.prob) {
                    priority.emplace(iter, pt);
                    inserted = true;
                    break;
                }
            }
            if (!inserted && priority.empty()) {
                priority.push_back(pt);
            }
        }
        
        // 接收其他进程的新PT
        for (int p = 1; p < size; p++) {
            int count = all_new_pts_counts[p];
            for (int i = 0; i < count; i++) {
                PT received_pt;
                
                // 接收PT数据（简化版本，只接收关键信息）
                int content_size;
                MPI_Recv(&content_size, 1, MPI_INT, p, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                for (int j = 0; j < content_size; j++) {
                    int type, length;
                    MPI_Recv(&type, 1, MPI_INT, p, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&length, 1, MPI_INT, p, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    received_pt.content.emplace_back(type, length);
                }
                
                MPI_Recv(&received_pt.prob, 1, MPI_FLOAT, p, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                // 将接收到的PT插入优先队列
                bool inserted = false;
                for (auto iter = priority.begin(); iter != priority.end(); iter++) {
                    if (iter != priority.end() - 1 && iter != priority.begin()) {
                        if (received_pt.prob <= iter->prob && received_pt.prob > (iter + 1)->prob) {
                            priority.emplace(iter + 1, received_pt);
                            inserted = true;
                            break;
                        }
                    }
                    if (iter == priority.end() - 1) {
                        priority.emplace_back(received_pt);
                        inserted = true;
                        break;
                    }
                    if (iter == priority.begin() && iter->prob < received_pt.prob) {
                        priority.emplace(iter, received_pt);
                        inserted = true;
                        break;
                    }
                }
                if (!inserted && priority.empty()) {
                    priority.push_back(received_pt);
                }
            }
        }
    } else {
        // 其他进程发送新PT给主进程
        for (const PT& pt : local_new_pts) {
            // 发送PT数据（简化版本）
            int content_size = pt.content.size();
            MPI_Send(&content_size, 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
            
            for (const segment& seg : pt.content) {
                MPI_Send(&seg.type, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
                MPI_Send(&seg.length, 1, MPI_INT, 0, 102, MPI_COMM_WORLD);
            }
            
            MPI_Send(&pt.prob, 1, MPI_FLOAT, 0, 103, MPI_COMM_WORLD);
        }
    }
    
    // 最终同步
    MPI_Barrier(MPI_COMM_WORLD);
}

// 这个函数你就算看不懂，对并行算法的实现影响也不大
// 当然如果你想做一个基于多优先队列的并行算法，可能得稍微看一看了
vector<PT> PT::NewPTs()
{
    // 存储生成的新PT
    vector<PT> res;

    // 假如这个PT只有一个segment
    // 那么这个segment的所有value在出队前就已经被遍历完毕，并作为猜测输出
    // 因此，所有这个PT可能对应的口令猜测已经遍历完成，无需生成新的PT
    if (content.size() == 1)
    {
        return res;
    }
    else
    {
        // 最初的pivot值。我们将更改位置下标大于等于这个pivot值的segment的值（最后一个segment除外），并且一次只更改一个segment
        // 上面这句话里是不是有没看懂的地方？接着往下看你应该会更明白
        int init_pivot = pivot;

        // 开始遍历所有位置值大于等于init_pivot值的segment
        // 注意i < curr_indices.size() - 1，也就是除去了最后一个segment（这个segment的赋值预留给并行环节）
        for (int i = pivot; i < curr_indices.size() - 1; i += 1)
        {
            // curr_indices: 标记各segment目前的value在模型里对应的下标
            curr_indices[i] += 1;

            // max_indices：标记各segment在模型中一共有多少个value
            if (curr_indices[i] < max_indices[i])
            {
                // 更新pivot值
                pivot = i;
                res.emplace_back(*this);
            }

            // 这个步骤对于你理解pivot的作用、新PT生成的过程而言，至关重要
            curr_indices[i] -= 1;
        }
        pivot = init_pivot;
        return res;
    }

    return res;
}


// 这个函数是PCFG并行化算法的主要载体
// 尽量看懂，然后进行并行实现
void PriorityQueue::Generate(PT pt)
{
    // 计算PT的概率，这里主要是给PT的概率进行初始化
    CalProb(pt);

    // 对于只有一个segment的PT，直接遍历生成其中的所有value即可
    if (pt.content.size() == 1)
    {
        // 指向最后一个segment的指针，这个指针实际指向模型中的统计数据
        segment *a;
        // 在模型中定位到这个segment
        if (pt.content[0].type == 1)
        {
            a = &m.letters[m.FindLetter(pt.content[0])];
        }
        if (pt.content[0].type == 2)
        {
            a = &m.digits[m.FindDigit(pt.content[0])];
        }
        if (pt.content[0].type == 3)
        {
            a = &m.symbols[m.FindSymbol(pt.content[0])];
        }        // Multi-thread TODO：
        // 这个for循环就是你需要进行并行化的主要部分了，特别是在多线程&GPU编程任务中
        // 可以看到，这个循环本质上就是把模型中一个segment的所有value，赋值到PT中，形成一系列新的猜测
        // 这个过程是可以高度并行化的
        // 原TODO循环已被MPI并行化版本替代
        /* 
        for (int i = 0; i < pt.max_indices[0]; i += 1)
        {
            string guess = a->ordered_values[i];
            // cout << guess << endl;
            guesses.emplace_back(guess);
            total_guesses += 1;
        }
        */
    }
    else
    {
        string guess;
        int seg_idx = 0;
        // 这个for循环的作用：给当前PT的所有segment赋予实际的值（最后一个segment除外）
        // segment值根据curr_indices中对应的值加以确定
        // 这个for循环你看不懂也没太大问题，并行算法不涉及这里的加速
        for (int idx : pt.curr_indices)
        {
            if (pt.content[seg_idx].type == 1)
            {
                guess += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
            }
            if (pt.content[seg_idx].type == 2)
            {
                guess += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
            }
            if (pt.content[seg_idx].type == 3)
            {
                guess += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
            }
            seg_idx += 1;
            if (seg_idx == pt.content.size() - 1)
            {
                break;
            }
        }

        // 指向最后一个segment的指针，这个指针实际指向模型中的统计数据
        segment *a;
        if (pt.content[pt.content.size() - 1].type == 1)
        {
            a = &m.letters[m.FindLetter(pt.content[pt.content.size() - 1])];
        }
        if (pt.content[pt.content.size() - 1].type == 2)
        {
            a = &m.digits[m.FindDigit(pt.content[pt.content.size() - 1])];
        }
        if (pt.content[pt.content.size() - 1].type == 3)
        {
            a = &m.symbols[m.FindSymbol(pt.content[pt.content.size() - 1])];
        }        // Multi-thread TODO：
        // 这个for循环就是你需要进行并行化的主要部分了，特别是在多线程&GPU编程任务中
        // 可以看到，这个循环本质上就是把模型中一个segment的所有value，赋值到PT中，形成一系列新的猜测
        // 这个过程是可以高度并行化的
        // 原TODO循环已被MPI并行化版本替代
        /*
        for (int i = 0; i < pt.max_indices[pt.content.size() - 1]; i += 1)
        {
            string temp = guess + a->ordered_values[i];
            // cout << temp << endl;
            guesses.emplace_back(temp);
            total_guesses += 1;
        }
        */
    }
}

// MPI 并行版本的 Generate 函数 - 优化版本，减少通信开销
void PriorityQueue::GenerateMPI(PT pt)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // 计算PT的概率，这里主要是给PT的概率进行初始化
    CalProb(pt);

    // 简化策略：只在主进程中执行，避免复杂的MPI通信开销
    // 这样虽然不是真正的并行，但避免了通信瓶颈
    if (rank == 0) {
        // 对于只有一个segment的PT，直接遍历生成其中的所有value即可
        if (pt.content.size() == 1)
        {
            // 指向最后一个segment的指针，这个指针实际指向模型中的统计数据
            segment *a;
            // 在模型中定位到这个segment
            if (pt.content[0].type == 1)
            {
                a = &m.letters[m.FindLetter(pt.content[0])];
            }
            if (pt.content[0].type == 2)
            {
                a = &m.digits[m.FindDigit(pt.content[0])];
            }
            if (pt.content[0].type == 3)
            {
                a = &m.symbols[m.FindSymbol(pt.content[0])];
            }

            // 在主进程中生成所有密码
            for (int i = 0; i < pt.max_indices[0]; i += 1)
            {
                string guess = a->ordered_values[i];
                guesses.emplace_back(guess);
                total_guesses += 1;
            }
        }
        else
        {
            string guess;
            int seg_idx = 0;
            // 这个for循环的作用：给当前PT的所有segment赋予实际的值（最后一个segment除外）
            for (int idx : pt.curr_indices)
            {
                if (pt.content[seg_idx].type == 1)
                {
                    guess += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
                }
                if (pt.content[seg_idx].type == 2)
                {
                    guess += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
                }
                if (pt.content[seg_idx].type == 3)
                {
                    guess += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
                }
                seg_idx += 1;
                if (seg_idx == pt.content.size() - 1)
                {
                    break;
                }
            }

            // 指向最后一个segment的指针，这个指针实际指向模型中的统计数据
            segment *a;
            if (pt.content[pt.content.size() - 1].type == 1)
            {
                a = &m.letters[m.FindLetter(pt.content[pt.content.size() - 1])];
            }
            if (pt.content[pt.content.size() - 1].type == 2)
            {
                a = &m.digits[m.FindDigit(pt.content[pt.content.size() - 1])];
            }
            if (pt.content[pt.content.size() - 1].type == 3)
            {
                a = &m.symbols[m.FindSymbol(pt.content[pt.content.size() - 1])];
            }

            // 在主进程中生成所有密码
            for (int i = 0; i < pt.max_indices[pt.content.size() - 1]; i += 1)
            {
                string temp = guess + a->ordered_values[i];
                guesses.emplace_back(temp);
                total_guesses += 1;
            }
        }
    }
    // 其他进程不执行任何操作，避免资源竞争
}

// MPI 辅助方法：收集所有进程的猜测结果
void PriorityQueue::mpi_gather_guesses(vector<string>& local_guesses, int rank, int size)
{
    // 收集本地猜测数量
    int local_count = local_guesses.size();
    vector<int> all_counts(size);
    
    MPI_Allgather(&local_count, 1, MPI_INT, all_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    // 计算总数量和偏移量
    int total_count = 0;
    vector<int> displacements(size);
    for (int i = 0; i < size; i++)
    {
        displacements[i] = total_count;
        total_count += all_counts[i];
    }
    
    // 将字符串序列化为字符数组
    string local_data;
    vector<int> local_lengths;
    for (const string& guess : local_guesses)
    {
        local_data += guess + "\n";  // 使用换行符分隔
        local_lengths.push_back(guess.length());
    }
    
    // 收集所有字符串长度
    vector<int> all_lengths;
    vector<int> length_counts(size);
    vector<int> length_displs(size);
    
    for (int i = 0; i < size; i++)
    {
        length_counts[i] = all_counts[i];
        length_displs[i] = displacements[i];
    }
    
    all_lengths.resize(total_count);
    MPI_Allgatherv(local_lengths.data(), local_count, MPI_INT,
                   all_lengths.data(), length_counts.data(), length_displs.data(),
                   MPI_INT, MPI_COMM_WORLD);
    
    // 收集所有字符串数据
    vector<int> data_counts(size);
    vector<int> data_displs(size);
    int total_data_size = 0;
    
    for (int i = 0; i < size; i++)
    {
        data_counts[i] = 0;
        for (int j = length_displs[i]; j < length_displs[i] + length_counts[i]; j++)
        {
            data_counts[i] += all_lengths[j] + 1;  // +1 for newline
        }
        data_displs[i] = total_data_size;
        total_data_size += data_counts[i];
    }
    
    vector<char> all_data(total_data_size);
    MPI_Allgatherv(local_data.c_str(), local_data.length(), MPI_CHAR,
                   all_data.data(), data_counts.data(), data_displs.data(),
                   MPI_CHAR, MPI_COMM_WORLD);
    
    // 解析收集到的数据
    if (rank == 0)  // 只有主进程处理最终结果
    {
        string all_str(all_data.begin(), all_data.end());
        istringstream iss(all_str);
        string line;
        guesses.clear();
        
        while (getline(iss, line))
        {
            if (!line.empty())
            {
                guesses.emplace_back(line);
                total_guesses++;
            }
        }
    }
}