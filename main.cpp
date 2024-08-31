#include "global.h"
#include "graph.h"
#include "bitmap.h"
// #include "bplustree.h"
#include <memkind.h>
#include "tools.h"
#include "pma.h"
#include <set>
#include "benchmark.h"
#include <random>
#include <limits.h>

// #include "avltree.h"
using namespace std;

#define GRAPH_LAYOUT_NAME "graph"
#define MAX_BUF_LEN 31

#define POOL_SIZE (int32_t)1024 * 1024 * 1024  // 1 GB
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

//#define size_v int

vector<string> data_path;
vector<string> task_path;
void set_path(vector<string>& data_path, vector<string>& task_path)  
{
	data_path.push_back("./data/demo.txt");  // 0: demo 
	data_path.push_back("./data/demo1.txt");  // 1: demo1  edge list 
	data_path.push_back("/home/tongfeng/data/livejournal/out.el");  // 2:livejournal
	data_path.push_back("/home/tongfeng/data/roadnet-ca/out.el");  // 3:roadnet-ca
	data_path.push_back("/home/tongfeng/data/friendster/out.el");  // 4:friendster
	data_path.push_back("/home/tongfeng/data/dimacs10-uk/out.el");  // 5:dimacs10-uk
	data_path.push_back("/home/tongfeng/data/orkut-links/out.el");  // 6:orkut-links
	data_path.push_back("/home/tongfeng/data/skitter-undirect/out.el");  // 7:skitter
	data_path.push_back("/home/tongfeng/data/dblp-undirect/out.el");  // 8:dblp
	data_path.push_back("/home/tongfeng/data/graph500-26.e");  // 9:graph500-26.e 
	data_path.push_back("/home/tongfeng/data/dota-league.e");  // 10:dota-league.e 

	data_path.push_back("/home/tongfeng/data/graph500-22.e");  // 11:graph500-22.e 
	data_path.push_back("/home/tongfeng/data/graph500-22.task");  // 12:graph500-22.task 
	data_path.push_back("/home/tongfeng/data/uniform-22.task");  // 13: uniform-22.task
	data_path.push_back("/home/tongfeng/data/uniform-22.e");  // 14: 
}

void generate_tasks(vector<Task> &t_vec)
{
    
	int32_t task_cnt = 0;
	for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
		Vertex *v = global_vectex_vec[vidx];
		if (v->degree > 0) {
			for (int blk_idx = 0; blk_idx < v->blk_list.size(); blk_idx++) {
                for (int ele_idx = 0; ele_idx < PM_BLK_SIZE; ele_idx++) {
                    if (v->blk_list[blk_idx].blk_ptr[ele_idx] >= 0) {
						// cout << "neighbor: " << v->blk_list[blk_idx].blk_ptr[ele_idx] << endl;
						Task task;
                        new_edge e;
						e.src = vidx;
						e.des = v->blk_list[blk_idx].blk_ptr[ele_idx];
						task.e = e;
						task.task_id = global_task_id;
						task.type = 1;  // 
						t_vec.push_back(task);
						task_cnt++;
						break;
                    } 
                }
				break;  
            }
		} else {
			continue;
		}

		if (task_cnt % 100 == 0) { 
			Task task;
			task.type = 2; // pagerank
			task.task_id = global_task_id;
			t_vec.push_back(task);
			global_task_id++;
			// task_cnt++;
		}
		if (task_cnt == 400) {
			break;
		}
	}
	// cout << "task_cnt: " << task_cnt << endl;
	task_cnt = t_vec.size();
	for (int32_t taskIdx = 0; taskIdx < task_cnt; taskIdx++) {
		if (t_vec[taskIdx].type == 1) {
			Task task = t_vec[taskIdx];
			task.type = 0;
			task.task_id = global_task_id;
			t_vec.push_back(task);
		} else {
			Task task = t_vec[taskIdx];
			task.type = 2;
			task.task_id = global_task_id;
			t_vec.push_back(task);
			global_task_id++;
		}
	}

	// print task info
	cout << "total number of tasks: " << t_vec.size() << endl;
	// for (int32_t tIdx = 0; tIdx < t_vec.size(); tIdx++) {
	// 	cout << t_vec[tIdx].type << ", ";
	// }
}

void test_concurrency_benchmark(vector<Task> t_vec, int32_t pool_threads, memkind_t kind)
{
	double start_t, end_t;
	vector<std::future<void>> futures;
	ThreadPool pool(pool_threads);
	start_t = get_current_time();
	new_edge e;
	// futures.push_back(pool.enqueue(clean_snap_cur));

	for (int tIdx = 0; tIdx < t_vec.size(); tIdx++) {
		Task task = t_vec[tIdx];
		// new_src_set.insert(e.src);
		if (task.type != 2) {
			if (global_task_map.find(task.task_id) != global_task_map.end()) {
				global_task_map[task.e.src].insert(task.e.src);
			} else {
				set<int32_t> tSet;
				tSet.insert(task.e.src);
				global_task_map[task.e.src] = tSet;
			}
		}

		// cout << "task idx: " << tIdx << " task.type: " << task.type << endl;
		while (writerCnt > 0 && task.type == 2); 
		switch (task.type) {
			case 0:
				// insert edge
				// new_edge e;
				// e.src = task.src;
				// e.des = task.des;
				// cout << "Insertion task!" << endl;
				futures.push_back(pool.enqueue(insert_edge_for_concurrent, task.e, kind));
				// task.task_id = global_task_id;
				// futures.push_back(pool.enqueue(PageRankPullGS_mv2, task));
				break;
			case 1:
				// cout << "Deletion task!" << endl;
				// delete edge
				// new_edge e;
				// e.src = task.src;
				// e.des = task.des;
				// cout << "deletion- src: " << task.e.src << " des: " << task.e.des << endl; 
				futures.push_back(pool.enqueue(delete_edge_for_concurrent, task.e, kind));
				break;
			case 2:
				// pagerank
				cout << "Query task!" << endl;
				futures.push_back(pool.enqueue(PageRankPullGS_mv2, task));
				break;
			default:
				break;
		}
	}
	futures.push_back(pool.enqueue(clean_snap_cur));

	cout << "All tasks are submitted." << endl;
	
    for(auto& future : futures) {
        future.wait();
    }

	end_t = get_current_time();
  	cout << "Concurrency benchmark time: " << end_t - start_t << endl;
}

void test_func(memkind_t kind)
{
	// pm_blk blk = generate_blk(kind);  
	pm_blk* blk = new pm_blk(generate_blk(kind));
	global_vectex_vec[1]->blk_list[0].next_blk = blk;
	blk->task_id = 3;
	cout << "v1 next blk task id 1.4: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
}

int main(int argc, char* argv[])
{
	// std::cout << "Size of int32_t: " << sizeof(int32_t) << " bytes. Max: " << INT_MAX << std::endl; // 2147483647
	// std::cout << "Size of int: " << sizeof(int) << " bytes" << std::endl;
	// std::cout << "Size of uint32_t: " << sizeof(uint32_t) << " bytes. Max: " << UINT_MAX << std::endl; // 4294967295

	// cout << "size of int*: " << sizeof(int32_t*) << " bytes." << endl;
	// cout << "size of next ptr: " << sizeof(PM_Block*) << " bytes." << endl;
	// return 1;

	double start_t, end_t;
	struct memkind *pmem_kind;
    int err;
    // Create PMEM memory pool with specific size

    err = memkind_create_pmem(PATH, 0, &pmem_kind);
	if (err)
		cout << "NVM memory initial error!" << endl;
	// pmem_kind = MEMKIND_DEFAULT;

	// int32_t *blk_pool = (int32_t *)memkind_malloc(pmem_kind, 10000 * sizeof(int32_t));
	// for (size_t i=0; i < 10000; i++) {
	// 	cout << "i: " << i << " address: " << &blk_pool[i] << endl;

	// 	if (memkind_detect_kind(&blk_pool[i]) == MEMKIND_DEFAULT) {
    //         cout << "is dram" << endl;
    //     }
	// }

	// pm_blk blk = generate_blk(pmem_kind);
	// if (blk.next_blk == NULL) 
	// 	cout << "next blk is NULL" << endl;
	// cout << "task id: " << blk.task_id << endl;


	// return 1;

	// const char* path = "/optane/tongfeng/test_pool";
	WorkerParams params; // input
	set_path(data_path, task_path);

	struct timeval start, end;
	struct timezone tzp;

	string dataPath = data_path[atoi(argv[1])];
    params.input_path = dataPath;

	load_graph_dram(params); 

	cout << "number of edges: " << global_edge_vec.size() << " global_max_id: " << global_max_id << endl;
	// cout << "global_edge_num: " << global_edge_num << endl;
	// print_edges();

	graph_init(global_max_id, pmem_kind);

	// while(1);
	// std::sort(global_edge_vec.begin(), global_edge_vec.end(), [](new_edge &a, new_edge &b) {
    //     return a.src < b.src;
	// });

	// print_edges();

	vector<int64_t> offset_vec;
	// offset_vec.push_back(0);  
	// int32_t cur_vid = global_edge_vec[0].src;

	global_edge_array = (new_edge *)memkind_malloc(pmem_kind, sizeof(new_edge) * global_edge_num);
	// #pragma omp parallel for
	// for (int64_t i = 0; i < global_edge_num; i++) {
	// 	global_edge_array[i] = global_edge_vec[i];
	// }

	// for (int64_t i = 0; i < global_edge_num; i++) {
	// 	// global_edge_array[i] = global_edge_vec[i];
	// 	if (global_edge_vec[i].src != cur_vid) {
	// 		offset_vec.push_back(i);
	// 		cur_vid = global_edge_vec[i].src;
	// 	}
	// }
	// global_edge_vec.clear();
	// global_edge_vec.shrink_to_fit();
	// cout << "move to NVM done." << endl;

	// for (int32_t i = 0; i < offset_vec.size(); i++) {
	// 	if (offset_vec[i] >= global_edge_num || offset_vec[i] < 0) {
	// 		cout << "index error!" << endl;
	// 	}
	// }

	// cout << offset_vec[0] << " " << offset_vec[1] << " " << offset_vec[2] << " " << offset_vec[3] << endl;

	// graph_maintenance_nvm_parallel_shuffle(pmem_kind);  
	graph_maintenance_nvm_parallel_shuffle_mix(pmem_kind);
	
	int64_t total_blk_cnt = 0;
	for (int32_t i = 0; i < global_vectex_vec.size(); i++) {
		total_blk_cnt += global_vectex_vec[i]->blk_list.size();
	}
	cout << "total number of blocks: " << total_blk_cnt << endl;
	// graph_maintenance_nvm(pmem_kind);  
	// graph_maintenance_delete(pmem_kind);
	// test_benchmark();

	
	// vector<Task> t_vec;
	// generate_tasks(t_vec);
	// test_concurrency_benchmark(t_vec, (32 - THREAD_NUM), pmem_kind);
	return 1;

	int32_t stride = global_edge_num;
	int32_t start_idx = 0;
	int32_t end_idx = start_idx + stride;
	// Bplus* bplus = new Bplus();
	cout << "Load graph on NVM parallel begin!" << endl;
	start_t = get_current_time();
	double total_time = 0.0;


	omp_set_num_threads(64);
	start_t = get_current_time();
	#pragma omp parallel for
	for (int64_t i = 0; i < offset_vec.size(); i++) {
		// if (i == (offset_vec.size()-2)) {
		// 	cout << "The last 2!!" << endl;
		// }
		int64_t start_off = offset_vec[i];
		int64_t end_off;
		if (i == (offset_vec.size()-1)) {
			end_off = global_edge_num;  
		} else {
			end_off = offset_vec[i+1];
		}
		Vertex *src = global_vectex_vec[global_edge_array[start_off].src];
		// graph_maintenance_nvm_insertion_parallel(src, start_off, end_off, pmem_kind);

		graph_maintenance_nvm_mixed_parallel(src, start_off, end_off, pmem_kind);
		
		// if (i % 1000000 == 0) {
		// 	cout << "insert " << i << " edges." << endl;
		// }
		// graph_maintenance_nvm_parallel(src, thread_id, pmem_kind);  
	}

	end_t = get_current_time();
	total_time += (end_t - start_t);

	// print_graph();


	// while (end_idx <= global_edge_num){
	// 	
	// 	set<int32_t> new_src_set;
	// 	vector<int32_t> new_src_vec;
	// 	for (size_t edgeIdx = start_idx; edgeIdx < end_idx; edgeIdx++) {
	// 		new_edge e = global_edge_array[edgeIdx];
	// 		// cout << "(" << e.src << ", " << e.des << ") " << global_vectex_vec[e.src]->id << endl;
	// 		global_vectex_vec[e.src]->edge_stamp_vec[0].push_back(e.des); 
	// 		new_src_set.insert(e.src);
	// 	}

	// 	// #pragma omp parallel for
	// 	// for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
	// 	// 	int32_t new_edge_num = global_vectex_vec[vidx]->edge_stamp_vec[0].size();
	// 	// 	global_vectex_vec[vidx]->new_edge_cnt = new_edge_num;
	// 	// 	global_vectex_vec[vidx]->new_edge_nvm = (int32_t *)memkind_malloc(pmem_kind, sizeof(int32_t) * new_edge_num);
	// 	// 	for (int32_t eIdx = 0; eIdx < new_edge_num; eIdx++) {
	// 	// 		global_vectex_vec[vidx]->new_edge_nvm[eIdx] = global_vectex_vec[vidx]->edge_stamp_vec[0][eIdx];
	// 	// 	}
	// 	// 	global_vectex_vec[vidx]->edge_stamp_vec[0].clear();
	// 	// 	global_vectex_vec[vidx]->edge_stamp_vec[0].shrink_to_fit();
	// 	// }


	// 	// #pragma omp parallel for
	// 	for (const int& element : new_src_set) { 
    //     	// Vertex *src = global_vectex_vec[element];
	// 		new_src_vec.push_back(element);
	// 		// graph_maintenance_nvm_parallel(src, pmem_kind);  
    // 	}

	// 	// cout << "start_idx 0: " << start_idx << " new_src_vec size: " << new_src_vec.size() << endl;
	// 	cout << "new_src_vec size: " << new_src_vec.size() << endl;
	// 	omp_set_num_threads(THREAD_NUM);
	// 	start_t = get_current_time();
	// 	#pragma omp parallel for
	// 	for (int i=0; i<new_src_vec.size(); i++) {
	// 		// int thread_id = omp_get_thread_num();
	// 		// cout << "thread_id: " << thread_id << endl;
	// 		Vertex *src = global_vectex_vec[new_src_vec[i]];
	// 		graph_maintenance_nvm_parallel_baseline(src, pmem_kind);  
	// 		// graph_maintenance_nvm_parallel_forLarge(src, pmem_kind);
			
			
	// 		// graph_maintenance_nvm_parallel(src, thread_id, pmem_kind);  
	// 	}

	// 	end_t = get_current_time();
	// 	total_time += (end_t - start_t);
		
	// 	start_idx = end_idx;
	// 	if (end_idx < global_edge_num)
	// 	{
	// 		end_idx = MIN(start_idx + stride, global_edge_num);
	// 	} else {
	// 		break;
	// 	}
	// 	// cout << "start_idx 1: " << start_idx << endl;
	// }
	// end_t = get_current_time();
    // print_graph();
	// test_priority_traversal();
    // print_pma(global_vectex_vec[1]);
	// cout << "Load graph on NVM parallel time: " << end_t - start_t << endl;
	// global_blk_vec.clear();
	// global_blk_vec.shrink_to_fit();  
	cout << "Load graph on NVM parallel time: " << total_time << endl;
	test_benchmark();

	
	// int64_t total_blk_cnt = 0;
	for (int32_t i = 0; i < global_vectex_vec.size(); i++) {
		total_blk_cnt += global_vectex_vec[i]->blk_list.size();
	}
	cout << "total number of blocks: " << total_blk_cnt << endl;

	// while (1);
	// free(bplus);
	// bplus->treeTraversal();

	// struct memkind *pmem_kind1;
    // int err1;
    // // Create PMEM memory pool with specific size
    // err1 = memkind_create_pmem(PATH, 0, &pmem_kind1);

	// if (err1)
	// 	cout << "NVM memory initial error!" << endl;


	return 1;


	cout << "Load graph on NVM finished!" << endl;
	
	return 1;

	// graph_blk* csr_blk = (graph_blk*)memkind_malloc(pmem_kind, sizeof(graph_blk));;

	// sorted_csr_pm(params, csr_blk);
	// return 1;

	// int32_t total_blks = 1000000;
	// generate_blk_batch(total_blks, pmem_kind);
	// int32_t *tmp_buf = (int32_t *)memkind_malloc(pmem_kind, PM_BLK_SIZE * sizeof(int32_t) * total_blks);
	// random_device rd;
	// mt19937 gen(rd());
	// shuffle(global_blk_vec.begin(), global_blk_vec.end(), gen); 
	// start_t = get_current_time();
	// for (int32_t i = 0; i < (total_blks * PM_BLK_SIZE); i++) {
	// 	int32_t a = tmp_buf[i];
	// }
	// end_t = get_current_time();
	// cout << "read time 1: " << end_t - start_t << endl;

	// start_t = get_current_time();
	// #pragma omp parallel for
	// for (int32_t i = 0; i < global_blk_vec.size(); i++) {
	// 	// __builtin_prefetch (global_blk_vec[i].blk_ptr, 0, 3);
	// 	#pragma omp parallel for
	// 	for (int j = 0; j < PM_BLK_SIZE; j++) {
	// 		int a = global_blk_vec[i].blk_ptr[j];
	// 	}
	// }
	// end_t = get_current_time();
	// cout << "read time 2: " << end_t - start_t << endl;

    // return 0;

	// std::cout << "Size of unsigned long: " << sizeof(unsigned long int) << " bytes" << std::endl;
    // std::cout << "Size of uint32_t: " << sizeof(uint32_t) << " bytes" << std::endl;
	// std::cout << "Size of int32_t: " << sizeof(int32_t) << " bytes" << std::endl;
	// return 1;

	// int *a = (int *)memkind_malloc(pmem_kind, 10 * sizeof(int));
	// memset(a, -1, 10 * sizeof(int));
	// // for (int i = 0; i < 4; i++) {
	// // 	cout << a[i] << ", ";
	// // }
	// int *b = &a[0];
	// int *c = &a[5];

	// cout << "c[3]: " << c[3] << " a[5]: " << a[5] << endl;

	// return 1;
}
