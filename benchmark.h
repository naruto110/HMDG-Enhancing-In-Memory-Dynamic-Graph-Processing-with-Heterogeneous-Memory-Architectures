#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include "global.h"
#include "pma.h"
#include "gapbs/pvector.h"
#include "gapbs/bitmap.h"
#include "gapbs/timer.h"
#include "gapbs/sliding_queue.h"
#include "gapbs/platform_atomics.h"
#include "gapbs/util.h"
#include <iostream>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <chrono>

using namespace std; 

typedef float ScoreT;
const float kDamp = 0.85;
ScoreT * PageRankPullGS(int max_iters, double epsilon = 1e-4) 
{
  omp_set_num_threads(THREAD_NUM);
  const ScoreT init_score = 1.0f / global_max_id;
  const ScoreT base_score = (1.0f - kDamp) / global_max_id;
  ScoreT *scores;
  ScoreT *outgoing_contrib;

  scores = (ScoreT *) malloc(sizeof(ScoreT) * global_max_id);
  outgoing_contrib = (ScoreT *) malloc(sizeof(ScoreT) * global_max_id);

  #pragma omp parallel for
  // #pragma omp parallel for num_threads(THREAD_NUM)
  for (int32_t vidx = 0; vidx < global_max_id; vidx++)
    outgoing_contrib[vidx] = init_score;


 for (int iter=0; iter < max_iters; iter++) {
    double error = 0;
    // #pragma omp parallel for
    // for (int32_t n=0; n < global_max_id; n++)
    //   outgoing_contrib[n] = scores[n] / global_vectex_vec[n]->degree;

    // int32_t total_num_access = 0;
    
    // #pragma omp parallel for reduction(+ : error) schedule(dynamic, THREAD_NUM)
    #pragma omp parallel for reduction(+ : error) schedule(dynamic, THREAD_NUM)
    for (int32_t vidx=0; vidx < global_max_id; vidx++) {
      // printf("I am thread %d / %d \n", omp_get_thread_num(),omp_get_num_threads());
      // cout << "thread: " << omp_get_thread_num() << " total threads: " << omp_get_num_threads();
      ScoreT incoming_total = 0;
      Vertex *v = global_vectex_vec[vidx];
      int32_t blk_cnt = v->blk_list.size();
      // #pragma omp parallel for reduction(+ : incoming_total) schedule(dynamic, 64)
      for (int32_t blk_idx = 0; blk_idx < blk_cnt; blk_idx++) {
        int32_t end_idx = PM_BLK_SIZE;
        if (blk_idx == blk_cnt - 1 && (v->pma->arr_cap % PM_BLK_SIZE != 0)) {
          end_idx = v->pma->arr_cap % PM_BLK_SIZE;
        } 
        int32_t *blk_ptr = v->blk_list[blk_idx].blk_ptr;
        // cout << "end_idx: " << end_idx << endl;
        // #pragma omp parallel for reduction(+ : incoming_total) schedule(dynamic, 64)
        for (int32_t ele_idx = 0; ele_idx < end_idx; ele_idx++) {
            // total_num_access++;
            if (blk_ptr[ele_idx] >= 0) {
                incoming_total += outgoing_contrib[blk_ptr[ele_idx]];
            }
        }
      }

      ScoreT old_score = scores[vidx];
      scores[vidx] = base_score + kDamp * incoming_total;
      error += fabs(scores[vidx] - old_score);
      outgoing_contrib[vidx] = scores[vidx] / v->degree;
    }
    // cout << "iter: " << iter << " total_num_access: " << total_num_access << endl;
  //  printf(" %2d    %lf\n", iter, error);
  //  if (error < epsilon)
  //    break;
  }

  return scores;
}

ScoreT * PageRankPullGS_mv(int max_iters, Task task, double epsilon = 1e-4) 
{
  omp_set_num_threads(THREAD_NUM);
  const ScoreT init_score = 1.0f / global_max_id;
  const ScoreT base_score = (1.0f - kDamp) / global_max_id;
  ScoreT *scores;
  ScoreT *outgoing_contrib;

  scores = (ScoreT *) malloc(sizeof(ScoreT) * global_max_id);
  outgoing_contrib = (ScoreT *) malloc(sizeof(ScoreT) * global_max_id);

  #pragma omp parallel for
  // #pragma omp parallel for num_threads(THREAD_NUM)
  for (int32_t vidx = 0; vidx < global_max_id; vidx++)
    outgoing_contrib[vidx] = init_score;


 for (int iter=0; iter < max_iters; iter++) {
    double error = 0;
    // #pragma omp parallel for
    // for (int32_t n=0; n < global_max_id; n++)
    //   outgoing_contrib[n] = scores[n] / global_vectex_vec[n]->degree;

    // int32_t total_num_access = 0;
    
    // #pragma omp parallel for reduction(+ : error) schedule(dynamic, THREAD_NUM)  // 16384
    #pragma omp parallel for reduction(+ : error) schedule(dynamic, THREAD_NUM)
    for (int32_t vidx=0; vidx < global_max_id; vidx++) {
      // printf("I am thread %d / %d \n", omp_get_thread_num(),omp_get_num_threads());
      // cout << "thread: " << omp_get_thread_num() << " total threads: " << omp_get_num_threads();
      ScoreT incoming_total = 0;
      Vertex *v = global_vectex_vec[vidx];
      std::shared_lock<std::shared_timed_mutex> lock(v->mutex);  // 

      int32_t blk_cnt = v->blk_list.size();
      for (int32_t blk_idx = 0; blk_idx < blk_cnt; blk_idx++) {
        pm_blk *cur_blk = &(v->blk_list[blk_idx]);
        int32_t *blk_ptr = NULL;
        while (cur_blk->next_blk != NULL)
        {
          if (cur_blk->task_id > task.task_id) {
            blk_ptr = cur_blk->blk_ptr;
            break;
          }
          cur_blk = cur_blk->next_blk;
        }
        if (cur_blk->task_id <= task.task_id) {
          blk_ptr = v->blk_list[blk_idx].blk_ptr;
        } else { 
          blk_ptr = cur_blk->blk_ptr;
        }
        
        int32_t end_idx = PM_BLK_SIZE;
        if (blk_idx == blk_cnt - 1 && (v->pma->arr_cap % PM_BLK_SIZE != 0)) {
          end_idx = v->pma->arr_cap % PM_BLK_SIZE;
        } 
        // int32_t *blk_ptr = v->blk_list[blk_idx].blk_ptr;
        // cout << "end_idx: " << end_idx << endl;
        for (int32_t ele_idx = 0; ele_idx < end_idx; ele_idx++) {
            // total_num_access++;
            if (blk_ptr[ele_idx] >= 0) {
                incoming_total += outgoing_contrib[blk_ptr[ele_idx]];
            }
        }
      }

      ScoreT old_score = scores[vidx];
      scores[vidx] = base_score + kDamp * incoming_total;
      error += fabs(scores[vidx] - old_score);
      outgoing_contrib[vidx] = scores[vidx] / v->degree;
    }
    // cout << "iter: " << iter << " total_num_access: " << total_num_access << endl;
  //  printf(" %2d    %lf\n", iter, error);
  //  if (error < epsilon)
  //    break;
  }

  return scores;
}

void PageRankPullGS_mv2(Task task) 
{
  // cout << "test pagerank.." << endl;
  omp_set_num_threads(THREAD_NUM);
  int max_iters = 20;
  double epsilon = 1e-4;
  const ScoreT init_score = 1.0f / global_max_id;
  const ScoreT base_score = (1.0f - kDamp) / global_max_id;
  ScoreT *scores;
  ScoreT *outgoing_contrib;

  scores = (ScoreT *) malloc(sizeof(ScoreT) * global_max_id);
  outgoing_contrib = (ScoreT *) malloc(sizeof(ScoreT) * global_max_id);

  #pragma omp parallel for
  // #pragma omp parallel for num_threads(THREAD_NUM)
  for (int32_t vidx = 0; vidx < global_max_id; vidx++)
    outgoing_contrib[vidx] = init_score;

//  cout << "test --0" << endl;
 for (int iter=0; iter < max_iters; iter++) {
    // cout << "pagerank iter: " << iter << endl;
    double error = 0;
    // #pragma omp parallel for
    // for (int32_t n=0; n < global_max_id; n++)
    //   outgoing_contrib[n] = scores[n] / global_vectex_vec[n]->degree;

    // int32_t total_num_access = 0;
    
    // #pragma omp parallel for reduction(+ : error) schedule(dynamic, THREAD_NUM)  // 16384
    #pragma omp parallel for reduction(+ : error) schedule(dynamic, THREAD_NUM)
    for (int32_t vidx=0; vidx < global_max_id; vidx++) {
      // printf("I am thread %d / %d \n", omp_get_thread_num(),omp_get_num_threads());
      // cout << "thread: " << omp_get_thread_num() << " total threads: " << omp_get_num_threads();
      ScoreT incoming_total = 0;
      Vertex *v = global_vectex_vec[vidx];
      // std::shared_lock<std::shared_timed_mutex> lock(v->mutex);  // 

      int32_t blk_cnt = v->blk_list.size();
       
      for (int32_t blk_idx = 0; blk_idx < blk_cnt; blk_idx++) {
        pm_blk *cur_blk = &(v->blk_list[blk_idx]);
        int32_t *blk_ptr = NULL;
        while (cur_blk->next_blk != NULL)
        {
          if (cur_blk->task_id > task.task_id) {
            blk_ptr = cur_blk->blk_ptr;
            break;
          }
          cur_blk = cur_blk->next_blk;
        }
        if (cur_blk->task_id <= task.task_id) {
          blk_ptr = v->blk_list[blk_idx].blk_ptr;
        } else { 
          blk_ptr = cur_blk->blk_ptr;
        }
        
        int32_t end_idx = PM_BLK_SIZE;
        if (blk_idx == blk_cnt - 1 && (v->pma->arr_cap % PM_BLK_SIZE != 0)) {
          end_idx = v->pma->arr_cap % PM_BLK_SIZE;
        } 
        // int32_t *blk_ptr = v->blk_list[blk_idx].blk_ptr;
        // cout << "end_idx: " << end_idx << endl;
        for (int32_t ele_idx = 0; ele_idx < end_idx; ele_idx++) {
            // total_num_access++;
            if (blk_ptr[ele_idx] >= 0) {
                incoming_total += outgoing_contrib[blk_ptr[ele_idx]];
            }
        }
      }

      ScoreT old_score = scores[vidx];
      scores[vidx] = base_score + kDamp * incoming_total;
      error += fabs(scores[vidx] - old_score);
      outgoing_contrib[vidx] = scores[vidx] / v->degree;
    }
    // cout << "iter: " << iter << " total_num_access: " << total_num_access << endl;
  //  printf(" %2d    %lf\n", iter, error);
  //  if (error < epsilon)
  //    break;
  }
  completeTaskVec[task.task_id] = 1;
}


void test_priority_traversal(void)
{
    double start_t = get_current_time();
    int32_t priority_edges = 0;
    #pragma omp parallel for schedule(dynamic, 1048576) 
    for (int32_t vidx = 0; vidx < global_max_id; vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        int32_t blk_num = v->blk_list.size();
        for (int blk_idx = (blk_num-1); blk_idx >= 0; blk_idx--) {
          int32_t end_idx = PM_BLK_SIZE;
          if (blk_idx == blk_num - 1 && (v->pma->arr_cap % PM_BLK_SIZE != 0)) {
            end_idx = v->pma->arr_cap % PM_BLK_SIZE;
          } 
            int32_t *blk_ptr = v->blk_list[blk_idx].blk_ptr;
            bool break_flag = false;
            for (int ele_idx = (end_idx-1); ele_idx >= 0; ele_idx--) {
                if (blk_ptr[ele_idx] > vidx) {
                    priority_edges++;
                } else if (blk_ptr[ele_idx] >= 0 && blk_ptr[ele_idx] < vidx){
                    break_flag = true;
                    break;
                }
            }
            if (break_flag)
                break;
        }
          // cout << endl;
      
    }
    cout << "total number of edges: " << priority_edges << endl;
    double end_t = get_current_time();
	cout << "Priority traversal time: " << end_t - start_t << endl;
}

// BFS

int32_t BUStep(pvector<int32_t> &parent, Bitmap &front,
               Bitmap &next) {
  int32_t awake_count = 0;
  next.reset();
  #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
  for (int32_t uidx=0; uidx < global_max_id; uidx++) {
    if (parent[uidx] < 0) {
      Vertex *u = global_vectex_vec[uidx];
      int32_t blk_num = u->blk_list.size();
      int32_t arr_cap = u->pma->arr_cap;
      for (int blk_idx = (blk_num-1); blk_idx >= 0; blk_idx--) {
        int32_t end_idx = PM_BLK_SIZE;
        if (blk_idx == blk_num - 1 && (arr_cap % PM_BLK_SIZE != 0)) {
          end_idx = arr_cap % PM_BLK_SIZE;
        } 
          int32_t *blk_ptr = u->blk_list[blk_idx].blk_ptr;
          for (int ele_idx = (end_idx-1); ele_idx >= 0; ele_idx--) {
            int32_t nei_v = blk_ptr[ele_idx];
              if (nei_v >= 0) {
                  if (front.get_bit(nei_v)) {
                    parent[uidx] = nei_v;
                    awake_count++;
                    next.set_bit(uidx);
                    break;
                  }
              }
          }
      }
    }
  }
  return awake_count;
}


int64_t TDStep(pvector<int32_t> &parent, SlidingQueue<int32_t> &queue) {
  int64_t scout_count = 0;
  #pragma omp parallel
  {
    QueueBuffer<int32_t> lqueue(queue);
    #pragma omp for reduction(+ : scout_count)
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
      int32_t uidx = *q_iter;
      Vertex *u = global_vectex_vec[uidx];
      int32_t blk_num = u->blk_list.size();
      int32_t arr_cap = u->pma->arr_cap;
      for (int blk_idx = (blk_num-1); blk_idx >= 0; blk_idx--) {
        int32_t end_idx = PM_BLK_SIZE;
        if (blk_idx == blk_num - 1 && (arr_cap % PM_BLK_SIZE != 0)) {
          end_idx = arr_cap % PM_BLK_SIZE;
        } 
        int32_t *blk_ptr = u->blk_list[blk_idx].blk_ptr;
        for (int ele_idx = (end_idx-1); ele_idx >= 0; ele_idx--) {
          int32_t nei_v = blk_ptr[ele_idx];
          if (nei_v >= 0) {
            int32_t curr_val = parent[nei_v];
            if (curr_val < 0) {
              if (compare_and_swap(parent[nei_v], curr_val, uidx)) {
                lqueue.push_back(nei_v);
                scout_count += -curr_val;
              }
            }
          }
        }
      }
    }
    lqueue.flush();
  }
  return scout_count;
}


void QueueToBitmap(const SlidingQueue<int32_t> &queue, Bitmap &bm) {
  #pragma omp parallel for
  for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
    int32_t u = *q_iter;
    bm.set_bit_atomic(u);
  }
}

void BitmapToQueue(const Bitmap &bm, SlidingQueue<int32_t> &queue) {
  #pragma omp parallel
  {
    QueueBuffer<int32_t> lqueue(queue);
    #pragma omp for
    for (int32_t n=0; n < global_max_id; n++)
      if (bm.get_bit(n))
        lqueue.push_back(n);
    lqueue.flush();
  }
  queue.slide_window();
}

pvector<int32_t> InitParent(void) {
  pvector<int32_t> parent(global_max_id);
  #pragma omp parallel for
  for (int32_t n=0; n < global_max_id; n++)
    parent[n] = global_vectex_vec[n]->degree != 0 ? -global_vectex_vec[n]->degree : -1;
  return parent;
}

pvector<int32_t> DOBFS(int32_t source, int alpha = 15, int beta = 18) {
//  PrintStep("Source", static_cast<int32_t>(source));
  Timer t;
  t.Start();
  pvector<int32_t> parent = InitParent();
  t.Stop();
//  PrintStep("i", t.Seconds());
  parent[source] = source;
  SlidingQueue<int32_t> queue(global_max_id);
  queue.push_back(source);
  queue.slide_window();
  Bitmap curr(global_max_id);
  curr.reset();
  Bitmap front(global_max_id);
  front.reset();
  int64_t edges_to_check = global_edge_num;
  int64_t scout_count = global_vectex_vec[source]->degree;
  while (!queue.empty()) {
    if (scout_count > edges_to_check / alpha) {
      int32_t awake_count, old_awake_count;
      TIME_OP(t, QueueToBitmap(queue, front));
//      PrintStep("e", t.Seconds());
      awake_count = queue.size();
      queue.slide_window();
      do {
        t.Start();
        old_awake_count = awake_count;
        awake_count = BUStep(parent, front, curr);
        front.swap(curr);
        t.Stop();
//        PrintStep("bu", t.Seconds(), awake_count);
      } while ((awake_count >= old_awake_count) ||
               (awake_count > global_max_id / beta));
      TIME_OP(t, BitmapToQueue(front, queue));
//      PrintStep("c", t.Seconds());
      scout_count = 1;
    } else {
      t.Start();
      edges_to_check -= scout_count;
      scout_count = TDStep(parent, queue);
      queue.slide_window();
      t.Stop();
//      PrintStep("td", t.Seconds(), queue.size());
    }
  }
  #pragma omp parallel for
  for (int32_t n = 0; n < global_max_id; n++)
    if (parent[n] < -1)
      parent[n] = -1;
  return parent;
}

// bc
static void ParallelPrefixSum(const pvector<int32_t> &degrees, pvector<int32_t> &prefix) {
  const size_t block_size = 1 << 20;
  const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
  pvector<int32_t> local_sums(num_blocks);
#pragma omp parallel for
  for (size_t block = 0; block < num_blocks; block++) {
    int32_t lsum = 0;
    size_t block_end = std::min((block + 1) * block_size, degrees.size());
    for (size_t i = block * block_size; i < block_end; i++) lsum += degrees[i];
    local_sums[block] = lsum;
  }
  pvector<int32_t> bulk_prefix(num_blocks + 1);
  int32_t total = 0;
  for (size_t block = 0; block < num_blocks; block++) {
    bulk_prefix[block] = total;
    total += local_sums[block];
  }
  bulk_prefix[num_blocks] = total;
#pragma omp parallel for
  for (size_t block = 0; block < num_blocks; block++) {
    int32_t local_total = bulk_prefix[block];
    size_t block_end = std::min((block + 1) * block_size, degrees.size());
    for (size_t i = block * block_size; i < block_end; i++) {
      prefix[i] = local_total;
      local_total += degrees[i];
    }
  }
  prefix[degrees.size()] = bulk_prefix[num_blocks];
}

static void ParallelLoadDegrees(pvector<int32_t> &degrees) {
#pragma omp parallel for schedule(dynamic, 16384)
  for (int32_t u = 0; u < global_max_id; u++) {
    degrees[u] = global_vectex_vec[u]->degree;
  }
}

inline int32_t GetEdgeId(const pvector<int32_t> &prefix, int32_t u, int32_t local_edge_id) {
  return prefix[u] + local_edge_id;
}

void PBFS(int32_t source, pvector<int32_t> &path_counts,
    Bitmap &succ, vector<SlidingQueue<int32_t>::iterator> &depth_index,
    SlidingQueue<int32_t> &queue, const pvector<int32_t> &prefix) {
  pvector<int32_t> depths(global_max_id, -1);
  depths[source] = 0;
  path_counts[source] = 1;
  queue.push_back(source);
  depth_index.push_back(queue.begin());
  queue.slide_window();

  #pragma omp parallel
  {
    int32_t depth = 0;
    QueueBuffer<int32_t> lqueue(queue);
    while (!queue.empty()) {
      #pragma omp single
      depth_index.push_back(queue.begin());
      depth++;
      #pragma omp for schedule(dynamic, 64)
      for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
        int32_t uidx = *q_iter;
        int32_t local_edge_id = 0;

        Vertex *u = global_vectex_vec[uidx];
        int32_t blk_num = u->blk_list.size();
        int32_t arr_cap = u->pma->arr_cap;
        for (int blk_idx = (blk_num-1); blk_idx >= 0; blk_idx--) {
          int32_t end_idx = PM_BLK_SIZE;
          if (blk_idx == blk_num - 1 && (arr_cap % PM_BLK_SIZE != 0)) {
            end_idx = arr_cap % PM_BLK_SIZE;
          } 
          int32_t *blk_ptr = u->blk_list[blk_idx].blk_ptr;
          for (int ele_idx = (end_idx-1); ele_idx >= 0; ele_idx--) {
            int32_t nei_v = blk_ptr[ele_idx];
            if (nei_v >= 0) {
              if ((depths[nei_v] == -1) && (compare_and_swap(depths[nei_v], static_cast<int32_t>(-1), depth))) {
                lqueue.push_back(nei_v);
              }
              if (depths[nei_v] == depth) {
                succ.set_bit_atomic(GetEdgeId(prefix, uidx, local_edge_id));
                fetch_and_add(path_counts[nei_v], path_counts[uidx]);
              }
              local_edge_id += 1;
            }
          }
        }
      }
      lqueue.flush();
      #pragma omp barrier
      #pragma omp single
      queue.slide_window();
    }
  }
  depth_index.push_back(queue.begin());
}


pvector<ScoreT> Brandes(int32_t source, int32_t num_iters) {
  pvector<ScoreT> scores(global_max_id, 0);
  pvector<int32_t> path_counts(global_max_id);
  Bitmap succ(global_edge_num);
  vector<SlidingQueue<int32_t>::iterator> depth_index;
  SlidingQueue<int32_t> queue(global_max_id);

  pvector<int32_t> degrees(global_max_id);
  ParallelLoadDegrees(degrees);
  pvector<int32_t> prefix(degrees.size() + 1);
  ParallelPrefixSum(degrees, prefix);

  for (int32_t iter=0; iter < num_iters; iter++) {
    int32_t source = source;
    path_counts.fill(0);
    depth_index.resize(0);
    queue.reset();
    succ.reset();
    PBFS(source, path_counts, succ, depth_index, queue, prefix);
    pvector<ScoreT> deltas(global_max_id, 0);
    for (int d=depth_index.size()-2; d >= 0; d--) {
      #pragma omp parallel for schedule(dynamic, 64)
      for (auto it = depth_index[d]; it < depth_index[d+1]; it++) {
        int32_t uidx = *it;
        ScoreT delta_u = 0;
        int32_t local_edge_id = 0;

        Vertex *u = global_vectex_vec[uidx];
        int32_t blk_num = u->blk_list.size();
        int32_t arr_cap = u->pma->arr_cap;
        for (int blk_idx = (blk_num-1); blk_idx >= 0; blk_idx--) {
          int32_t end_idx = PM_BLK_SIZE;
          if (blk_idx == blk_num - 1 && (arr_cap % PM_BLK_SIZE != 0)) {
            end_idx = arr_cap % PM_BLK_SIZE;
          } 
          int32_t *blk_ptr = u->blk_list[blk_idx].blk_ptr;
          for (int ele_idx = (end_idx-1); ele_idx >= 0; ele_idx--) {
            int32_t nei_v = blk_ptr[ele_idx];
            if (nei_v >= 0) {
              if (succ.get_bit(GetEdgeId(prefix, uidx, local_edge_id))) {
                delta_u += static_cast<ScoreT>(path_counts[uidx]) / static_cast<ScoreT>(path_counts[nei_v]) * (1 + deltas[nei_v]);
              }
              local_edge_id += 1;
            }
          }
        }

        deltas[uidx] = delta_u;
        scores[uidx] += delta_u;
      }
    }
  }
  // normalize scores
  ScoreT biggest_score = 0;
  #pragma omp parallel for reduction(max : biggest_score)
  for (int32_t n=0; n < global_max_id; n++)
    biggest_score = max(biggest_score, scores[n]);
  #pragma omp parallel for
  for (int32_t n=0; n < global_max_id; n++)
    scores[n] = scores[n] / biggest_score;
  return scores;
}


// cc
pvector<int32_t> ShiloachVishkin(void) {
  pvector<int32_t> comp(global_max_id);
  #pragma omp parallel for
  for (int32_t n=0; n < global_max_id; n++)
    comp[n] = n;
  bool change = true;
  int num_iter = 0;
  while (change) {
    change = false;
    num_iter++;
    // note: this gives better scaleup performance
    // #pragma omp parallel for schedule(dynamic, 64)
    #pragma omp parallel for
    for (int32_t uidx=0; uidx < global_max_id; uidx++) {
      Vertex *u = global_vectex_vec[uidx];
      int32_t blk_num = u->blk_list.size();
      int32_t arr_cap = u->pma->arr_cap;
      for (int blk_idx = (blk_num-1); blk_idx >= 0; blk_idx--) {
        int32_t end_idx = PM_BLK_SIZE;
        if (blk_idx == blk_num - 1 && (arr_cap % PM_BLK_SIZE != 0)) {
          end_idx = arr_cap % PM_BLK_SIZE;
        } 
        int32_t *blk_ptr = u->blk_list[blk_idx].blk_ptr;
        for (int ele_idx = (end_idx-1); ele_idx >= 0; ele_idx--) {
          int32_t nei_v = blk_ptr[ele_idx];
          if (nei_v >= 0) {
            int32_t comp_u = comp[uidx];
            int32_t comp_v = comp[nei_v];
            if (comp_u == comp_v) continue;
            // Hooking condition so lower component ID wins independent of direction
            int32_t high_comp = comp_u > comp_v ? comp_u : comp_v;
            int32_t low_comp = comp_u + (comp_v - high_comp);
            if (high_comp == comp[high_comp]) {
              change = true;
              comp[high_comp] = low_comp;
            }
          }
        }
      }
    }
    #pragma omp parallel for
    for (int32_t n=0; n < global_max_id; n++) {
      while (comp[n] != comp[comp[n]]) {
        comp[n] = comp[comp[n]];
      }
    }
  }
  cout << "Shiloach-Vishkin took " << num_iter << " iterations" << endl;
  return comp;
}

void test_benchmark(void)
{
    double start_t, end_t;
  	cout << "PageRank start." << endl;
		start_t = get_current_time();
		PageRankPullGS(10);
		end_t = get_current_time();
		cout << "PageRank time: " << end_t - start_t << endl;

    cout << "BFS start." << endl;
		start_t = get_current_time();
		DOBFS(1);
		end_t = get_current_time();
		cout << "BFS time: " << end_t - start_t << endl;

		// cout << "BC start." << endl;
		// start_t = get_current_time();
		// Brandes(1, 1);
		// end_t = get_current_time();
		// cout << "BC time: " << end_t - start_t << endl;
		
    // cout << "CC start." << endl;
		// start_t = get_current_time();
		// ShiloachVishkin();
		// end_t = get_current_time();
		// cout << "CC time: " << end_t - start_t << endl;
}

void test_PR_parallel(Task task)
{
  double start_t, end_t;

  // start_t = get_current_time();
  // #pragma omp parallel for
  // for (int i = 0; i < 3; i++) {
  //   PageRankPullGS_mv(20, task);
  // }
  // end_t = get_current_time();
  // cout << "OpenMP Parallel PageRank time: " << end_t - start_t << endl;

  cout << "PageRank start." << endl;
  start_t = get_current_time();
  PageRankPullGS(20);
  end_t = get_current_time();
  cout << "PageRank time: " << end_t - start_t << endl;

  cout << "PageRank start." << endl;
  start_t = get_current_time();
  PageRankPullGS_mv(20, task);
  end_t = get_current_time();
  cout << "PageRank_mv time: " << end_t - start_t << endl;

  cout << "PageRank start." << endl;
  start_t = get_current_time();
  // DOBFS(1);
  // DOBFS(1);
  // DOBFS(1);
  PageRankPullGS_mv(20, task);
  PageRankPullGS_mv(20, task);
  PageRankPullGS_mv(20, task);
  PageRankPullGS_mv(20, task);
  end_t = get_current_time();
  cout << "PageRank_mv time X3: " << end_t - start_t << endl;



  std::thread pr1([&](){
    PageRankPullGS_mv(20, task);
    // PageRankPullGS(20);
  });

  std::thread pr2([&](){
    PageRankPullGS_mv(20, task);
    // PageRankPullGS(20);
  });

  std::thread pr3([&](){
    PageRankPullGS_mv(20, task);
    // PageRankPullGS(20);
  });

    std::thread pr4([&](){
    PageRankPullGS_mv(20, task);
    // PageRankPullGS(20);
  });

  // start_t = get_current_time();
  // #pragma omp parallel
  // {
  //   #pragma omp single
  //   {
  //     #pragma omp task
  //     {
  //       PageRankPullGS_mv(20, task);
  //     }
  //     #pragma omp task
  //     {
  //       PageRankPullGS_mv(20, task);
  //     }
  //     #pragma omp task
  //     {
  //       PageRankPullGS_mv(20, task);
  //     }
  //   }
  // }
  // end_t = get_current_time();
  // cout << "OMP Parallel PageRank time: " << end_t - start_t << endl;

  start_t = get_current_time();
  pr1.join();
  pr2.join();
  pr3.join();
  pr4.join();
  end_t = get_current_time();
  cout << "Parallel PageRank time: " << end_t - start_t << endl;


}

class ThreadPool {
public:
    ThreadPool(size_t threads) : stop(false) {
        for(size_t i = 0; i<threads; ++i) {
            workers.emplace_back(
                [this] {
                    while(true) {
                        std::function<void()> task;

                        {
                            std::unique_lock<std::mutex> lock(this->queue_mutex);
                            this->condition.wait(lock, 
                                [this]{ return this->stop || !this->tasks.empty(); });
                            if(this->stop && this->tasks.empty())
                                return;
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                        }

                        task();
                    }
                }
            );
        }
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> std::future<typename std::result_of<F(Args...)>::type> {
        using return_type = typename std::result_of<F(Args...)>::type;

        auto task = std::make_shared< std::packaged_task<return_type()> >(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);

            // 不允许在停止的线程池中加入新任务
            if(stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            tasks.emplace([task](){ (*task)(); });
        }
        condition.notify_one();
        return res;
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for(std::thread &worker: workers)
            worker.join();
    }
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

void test_PR_threadpool(Task task)
{
  double start_t, end_t;
  cout << "PageRank start." << endl;
  start_t = get_current_time();
  PageRankPullGS_mv2(task);
  end_t = get_current_time();
  cout << "PageRank time: " << end_t - start_t << endl;

  ThreadPool pool(4);
  
  start_t = get_current_time();
  auto result1 = pool.enqueue(PageRankPullGS_mv2, task); 
  auto result2 = pool.enqueue(PageRankPullGS_mv2, task);

  std::cout << "Main thread continues executing other tasks..." << std::endl;

  auto result3 = pool.enqueue(PageRankPullGS_mv2, task);
  auto result4 = pool.enqueue(PageRankPullGS_mv2, task);

  result1.get();
  result2.get();
  result3.get();
  result4.get();

  end_t = get_current_time();
  cout << "PageRank time: " << end_t - start_t << endl;
}

#endif 
