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

#define size_v int32_t  
//#define MAX_ID 31+1  
#define BLK_SIZE 125  
#define MAX_LEN  1000000000  
#define MAX_DEG_CAP 20  
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
	int32_t *blk_ptr; 
	PM_Block* next_blk;
	// Bitmap32 myBitmap; 
} pm_blk;



typedef struct _pma
{
    int32_t ele_num;  
    int32_t arr_cap;  
    int32_t seg_len;  
    int32_t seg_num;  
    int32_t tree_h;  
    double delta_t;  
    double delta_p;  
} PMA;


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
	size_v next_blk;  
	size_v idle_slot;  
	size_v next_pos;  
	size_v buf[MAX_LEN];   // sizeof(block) = 4 KB  / 997
} graph_blk;


typedef struct one_block
{
	size_v* offset;  
	size_v* neiList;   
} one_blk;

struct WorkerParams {
	string input_path;
	string task_path;
	int global_mode_ctr;  // 0: PM only, 1: DRAM only, 2: PM+DRAM
};

typedef struct idleBlk
{
	size_v blk_idx;  
	size_v max_cap;  
	size_v next_pos;  
}idle_blk;

vector<graph_blk *> graph_blk_vec;  
vector<idle_blk *> idle_blk_vec;  // used in load graph
vector<idle_blk*> update_blk_vec;  // used in graph update
//Map idle_blk_map;  // key: blk_id, value: capacity
//MapIter mapIter;

graph_blk* generate_new_blk(PMEMoid root)
{
	graph_blk* blk = (graph_blk*)pmemobj_direct(root) + global_offset;  // global_offset * sizeof(graph_blk)
	global_offset += 1;
	blk->next_blk = -1; 
	blk->vcnt = 0;
	blk->blk_id = graph_blk_vec.size();
	blk->idle_slot = 0;
	blk->next_pos = MAX_LEN - 1;
	graph_blk_vec.push_back(blk);

	//cout << "blk id: " << blk->blk_id << endl;
	return blk;
}

int32_t global_max_id = 0;
int64_t global_edge_num = 0;
typedef pair<size_v, size_v> VIDX_PAIR;  // pair(blk_id, blk_offset)

// graph update
typedef struct new_edge
{
	int32_t src;  
	int32_t des;  
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
	bool is_process = 0; 
	int32_t edge_num = 0; 
	// int32_t new_nei_num = 0;
} remainSrc;

class Vertex {
public:
	int32_t id;
	int32_t degree;
	mutable std::shared_timed_mutex mutex;
	// int32_t blk_cnt;
	vector<pm_blk> blk_list;
	// pm_blk *blk_list; 
	PMA* pma;
	vector<vector<int32_t> > edge_stamp_vec;  
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
	int32_t type;  // 0: insert edge 1: delete edge 2: pagerank
} Task;


vector<new_edge> global_edge_vec;
new_edge *global_edge_array;
// int64_t global_new_edge_num; 
int32_t global_reserved_slots;

vector<pm_blk> global_blk_vec;

int32_t global_blk_id = 0;
vector<Vertex *> global_vectex_vec;
vector<int32_t *> global_array_vec;

double global_time_cnt = 0;
int32_t global_task_id = 0;
int32_t global_min_task = 0;  // for clean snapshot

#define THREAD_NUM 64
vector<vector<pm_blk> > global_parallel_blk_vec(THREAD_NUM);

// task id -> vid hashmap
taskMap global_task_map;

std::atomic<int> writerCnt(0);

std::array<int32_t, 100> completeTaskVec = {0};
// std::atomic<std::array<int32_t, 100>> completeTaskVec(initialData); 
// global_edge_record;

#endif

