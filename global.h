#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <libpmemobj.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include <ctime>
#include "fstream"
#include <sstream>
// #include <iostream>
#include <array>
#include "tools.h"
#include <queue>
#include "bplustree.h"
#include <omp.h>
// #include <algorithm>
// #include "pma.h"
#include <ext/hash_map>
#define hash_map __gnu_cxx::hash_map
#define hash_set __gnu_cxx::hash_set
#include "bitmap.h"

#include <shared_mutex>
#include <mutex>
#include <thread>
#include <atomic>
#include <set>

using namespace std;

#define size_v int32_t  // 4 bytes
//#define MAX_ID 31+1  // 最大id为31，+1是因为包含id为0的顶点
#define BLK_SIZE 125  // 4kb  int占4个字节 4000
#define MAX_LEN  1000000000  // 10亿slots
#define MAX_DEG_CAP 20  // 单个blk能装下最大degree顶点 4kb -> 826
#define BUF_PRO 0
#define PATH "/dcpmm/tongfeng"

size_v global_offset = 0;



typedef hash_map<int32_t, set<int32_t> > taskMap;
typedef typename taskMap::iterator taskIter;

/* Design for Dynamic Graph Processing on Heterogeneous Memory */
typedef struct PM_Block
{
	// int32_t blk_id;  
	int32_t task_id;
	// int32_t max_nei;
	// int32_t min_nei;
	int32_t *blk_ptr; // blk起始地址，nvm上的地址
	PM_Block* next_blk;
	// Bitmap32 myBitmap; // 记录当前blk各位置元素有效性，需要被初始化
} pm_blk;



typedef struct _pma
{
    int32_t ele_num;  // pma_array中元素个数
    int32_t arr_cap;  // pma_array size
    int32_t seg_len;  // segment length
    int32_t seg_num;  // segment 数量
    int32_t tree_h;   // tree height
    double delta_t;  // Delta for upper density threshold.
    double delta_p;  // Delta for lower density threshold
} PMA;

// typedef struct Vertex
// {
// 	int32_t id;
// 	int32_t degree;
// 	int32_t blk_cnt;
// 	pm_blk *blk_list; // 存储当前顶点邻节点blk链
// 	PMA* pma;
// } vertex;

/***************************************************************/

void* pm_alloc(const char *layout, int32_t poolsize)  // unit: byte
{
	const char* path = "/optane/tongfeng/test_pool";
	PMEMobjpool* pop = pmemobj_create(path, layout, poolsize, 0666); 
	PMEMoid root = pmemobj_root(pop, poolsize);
	void* mem_all = pmemobj_direct(root);
	return mem_all;
}

typedef struct block
{
	size_v blk_id;
	size_v vcnt;  // the number of vertices stored in the block
	size_v next_blk;  // 一个blk存不下，需要存储在另外的blk上
	size_v idle_slot;  // 当前blk还剩余多少空余slot
	size_v next_pos;  // 最后一个顶点最后一个slot后一个位置，用于计算最后一个顶点的buffer
	size_v buf[MAX_LEN];   // sizeof(block) = 4 KB  / 997
} graph_blk;


typedef struct one_block
{
	size_v* offset;  // 偏移量
	size_v* neiList;   // 邻节点数组，动态分配所需最大内存
} one_blk;

struct WorkerParams {
	string input_path;
	string task_path;
	int global_mode_ctr;  // 0: PM only, 1: DRAM only, 2: PM+DRAM
};

typedef struct idleBlk
{
	size_v blk_idx;  // blk在graph_blk_vec中的idex下标
	size_v max_cap;  // 最大degree的顶点(剩余slot数量)
	size_v next_pos;  //  gblk->buf最后一个位置的下一个位置
}idle_blk;

vector<graph_blk *> graph_blk_vec;  // 存储所有blk的首地址 block*
vector<idle_blk *> idle_blk_vec;  // 未存满的blk信息, used in load graph
vector<idle_blk*> update_blk_vec;  // used in graph update
//Map idle_blk_map;  // key: blk_id, value: capacity
//MapIter mapIter;

graph_blk* generate_new_blk(PMEMoid root)
{
	graph_blk* blk = (graph_blk*)pmemobj_direct(root) + global_offset;  // 注意：此处offset为指针offset，对应物理地址偏移量为 global_offset * sizeof(graph_blk)
	global_offset += 1;
	blk->next_blk = -1; // 初始化为-1，表明默认没有后续链接blk
	blk->vcnt = 0;
	blk->blk_id = graph_blk_vec.size();
	blk->idle_slot = 0;
	blk->next_pos = MAX_LEN - 1;
	graph_blk_vec.push_back(blk);

	//cout << "blk id: " << blk->blk_id << endl;
	return blk;
}

// 构建顶点索引
int32_t global_max_id = 0;
int64_t global_edge_num = 0;
typedef pair<size_v, size_v> VIDX_PAIR;  // pair(blk_id, blk_offset),根据顶点id索引对应所存储在的blk_id及偏移量

// graph update
typedef struct new_edge
{
	int32_t src;  // blk在graph_blk_vec中的idex下标
	int32_t des;  // 最大degree的顶点
	int statue; // 0: insertion <0 deletion
	// char type; // edge insertion, graph processing task...
} newEdge;

// graph update
typedef struct batch_edge
{
	size_v src;
	size_v edge_num; 
	bool is_process = false;
	Leaf_Node* l = NULL;
	// VAL_TYPE nei_avl = NULL;
} batchEdge;

typedef struct v_offset
{
	size_v index;
	size_v degree = 0; 
	VAL_TYPE nei_avl;
	int32_t *in_place_buf;  // 
	// int32_t new_nei_num = 0;
} vertexOff;

typedef struct remain_src
{
	size_v src_idx;
	int32_t src;
	bool is_process = 0; // 为0表示不需要进行window计算
	int32_t edge_num = 0; // 插入完成后置0
	// int32_t new_nei_num = 0;
} remainSrc;

class Vertex {
public:
	int32_t id;
	int32_t degree;
	mutable std::shared_timed_mutex mutex;
	// int32_t blk_cnt;
	vector<pm_blk> blk_list;
	// pm_blk *blk_list; // 存储当前顶点邻节点blk链
	PMA* pma;
	vector<vector<int32_t> > edge_stamp_vec;  // 统一批次中，按任务批次加入更新边，越后对应任务越新
	int32_t* new_edge_nvm;
	int32_t new_edge_cnt;

	Vertex(int32_t vid){
		id = vid;
		degree = 0;
		// blk_cnt = 1;
		// vertexesV.resize(1);
	}

	~Vertex(){}
};

typedef struct _task
{
	int32_t task_id;
	int32_t src;
	int32_t des;
	new_edge e;
	int32_t type;  // 0: insert edge 1: delete edge 2: pagerank (实验以pagerank为例)
} Task;

// 图数据边集合
vector<new_edge> global_edge_vec;
new_edge *global_edge_array;
// int64_t global_new_edge_num; // 更新边数量
int32_t global_reserved_slots;

vector<pm_blk> global_blk_vec;

int32_t global_blk_id = 0;
vector<Vertex *> global_vectex_vec;
vector<int32_t *> global_array_vec;

double global_time_cnt = 0;
int32_t global_task_id = 0;
int32_t global_min_task = 0;  // for clean snapshot，应该清理的task

#define THREAD_NUM 64
vector<vector<pm_blk> > global_parallel_blk_vec(THREAD_NUM);

// 构建task id -> vid hashmap
taskMap global_task_map;

std::atomic<int> writerCnt(0);

std::array<int32_t, 100> completeTaskVec = {0};
// std::atomic<std::array<int32_t, 100>> completeTaskVec(initialData); // 检测对应task是否被处理
// global_edge_record;

#endif

