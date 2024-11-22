#ifndef __PMA_H_
#define __PMA_H_

#include<iostream>
#include "global.h"
#include <memkind.h>
#include "tools.h"
#include <cassert>
#include "bitmap.h"

#include <ext/hash_map>
#define hash_map __gnu_cxx::hash_map
#define hash_set __gnu_cxx::hash_set

typedef hash_map<int32_t, vector<int32_t> > Map;
typedef typename Map::iterator MapIter;

using namespace std;
 
/* Height-based (as opposed to depth-based) thresholds. */
/* Upper density thresholds. */
static const double t_h = 0.75;  /* root. */
static const double t_0 = 1.00;  /* leaves. */
/* Lower density thresholds. */
static const double p_h = 0.50;  /* root. */
static const double p_0 = 0.25;  /* leaves. */

static const uint8_t max_sparseness =  2; // 1 / p_0;

#define IS_EMPTY_SLOT 0
#define IS_SINGLE_ELEMENT 1
#define IS_CONTINUE_ELEMENT 2
#define IS_RESERVED_SLOT 3
#define IS_VERTEX_ID 4

#define BLK_SCALE 4
#define PM_BLK_SIZE (64 * BLK_SCALE)  

pm_blk generate_blk(memkind_t kind);
void generate_blk_batch(int32_t blk_num, memkind_t kind);
PMA* pma_creat (memkind_t kind);
int32_t blk_search(pm_blk blk, int32_t des);
int32_t pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from_blk, int32_t to_blk, int32_t &index);
int32_t pma_insert(PMA* pma, Vertex *src, int32_t des, int32_t index, memkind_t kind);
int32_t pma_insert_mv(PMA* pma, Vertex *src, int32_t des, int32_t index, memkind_t kind);
void rebalance(PMA* pma, Vertex *src, int32_t index, memkind_t kind);
void rebalance_mv(PMA* pma, Vertex *src, int32_t index, memkind_t kind);
void spread (PMA* pma, Vertex *src, int32_t from, int32_t to, int32_t occupancy, memkind_t kind);
void resize (PMA* pma, Vertex *src, memkind_t kind);
static void compute_capacity (PMA* pma);
void print_vertex(Vertex *v);
void print_pma(Vertex *v);
void print_graph(void);
void print_vertex(Vertex *v);
void print_vertex_mv(int32_t vid, Task task);
void graph_maintenance_nvm(memkind_t kind);
void graph_maintenance_nvm_parallel(Vertex *src, int thread_num, memkind_t kind);
void resize_parallel(PMA* pma, Vertex *src, int thread_num, memkind_t kind);
void rebalance_parallel(PMA* pma, Vertex *src, int32_t index, int thread_num, memkind_t kind);
int32_t pma_insert_parallel(PMA* pma, Vertex *src, int32_t des, int32_t index, int thread_num, memkind_t kind);
void graph_maintenance_delete(memkind_t kind);
void insert_edge_for_concurrent(new_edge e, memkind_t kind);
void delete_edge_for_concurrent(new_edge e, memkind_t kind);
void recycle_blk_for_del(vector<Vertex *> vset);
void clean_snapshot(int32_t task_id, vector<Vertex *> vset);
void clean_snap_cur(void);

void graph_maintenance_nvm_parallel_shuffle(memkind_t kind);
void graph_maintenance_nvm_parallel_shuffle_mix(memkind_t kind);
int32_t binary_pma_find_mix(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to);

void graph_maintenance_nvm_mixed_parallel(Vertex *src, int64_t start_off, int64_t end_off, memkind_t kind);

pm_blk generate_blk(memkind_t kind)
{
    pm_blk blk;
    // blk.blk_id = global_blk_id;
    blk.blk_ptr = (int32_t *)memkind_malloc(kind, PM_BLK_SIZE * sizeof(int32_t)); 
    blk.next_blk = NULL;
    blk.task_id = 0;
    memset(blk.blk_ptr, -1, PM_BLK_SIZE * sizeof(int32_t));
    // global_blk_vec.push_back(blk);
    // global_blk_id++;
    return blk;
}

// void generate_blk_batch(int32_t blk_num, memkind_t kind)
// {
//     // cout << "generate blk batch!" << endl;
//     static int malloc_cnt = 0;
//     malloc_cnt++;
//     cout << "malloc_cnt: " << malloc_cnt << endl;
//     // blk_num = 10000;
//     blk_num = global_edge_vec.size() / PM_BLK_SIZE * 4;
//     int32_t *blk_pool = (int32_t *)memkind_malloc(kind, (PM_BLK_SIZE * sizeof(int32_t) * blk_num));

//     // memset(blk_pool, -1, PM_BLK_SIZE * sizeof(int32_t) * blk_num);
//     global_array_vec.push_back(blk_pool);
//     // cout << "blk_pool[0]: " << blk_pool[1000] << endl;
//     int32_t blk_start_off = 0;
//     for (int i = 0; i < blk_num; i++) {
//         pm_blk blk;
//         blk.blk_ptr = &blk_pool[blk_start_off]; 
//         global_blk_vec.push_back(blk);
//         cout << "blk_start_off: " << blk_start_off << " total size: " << (blk_num*PM_BLK_SIZE) << " malloc_cnt: " << malloc_cnt << endl;

//         // assert (blk_start_off <= blk_num*PM_BLK_SIZE);
//         if (memkind_detect_kind(&blk_pool[blk_start_off]) == MEMKIND_DEFAULT) {
//             cout << "is dram" << endl;
//         }
//         // cout << "generate test. blk[0]: " << blk.blk_ptr[0] << endl;
//         // initializeBitmap(&blk.myBitmap);
//         blk_start_off += PM_BLK_SIZE;
//     }
//     while(1);
// }

void generate_blk_batch(int32_t blk_num, memkind_t kind)
{
    // static int64_t blk_cnt = 0;
    // blk_cnt += blk_num;
    // cout << "generate blk batch!" << endl;
    // blk_num = 10000;
    // blk_num = global_edge_vec.size() / PM_BLK_SIZE * 4;
    for (int i=0; i < blk_num; i++) {
        pm_blk blk;
        blk.next_blk = NULL;
        blk.task_id = 0;
        blk.blk_ptr = (int32_t *)memkind_malloc(kind, (PM_BLK_SIZE * sizeof(int32_t)));
        memset(blk.blk_ptr, -1, PM_BLK_SIZE * sizeof(int32_t));
        global_blk_vec.push_back(blk);

        // for (size_t j=0; j < PM_BLK_SIZE; j++) {
        //     if (memkind_detect_kind(&blk.blk_ptr[j]) == MEMKIND_DEFAULT) {
        //         // cout << "is dram" << endl;
        //     }
	    // }
    }


    // int32_t *blk_pool = (int32_t *)memkind_malloc(kind, (PM_BLK_SIZE * sizeof(int32_t) * blk_num));

    // // memset(blk_pool, -1, PM_BLK_SIZE * sizeof(int32_t) * blk_num);
    // global_array_vec.push_back(blk_pool);
    // // cout << "blk_pool[0]: " << blk_pool[1000] << endl;
    // int32_t blk_start_off = 0;
    // for (int i = 0; i < blk_num; i++) {
    //     pm_blk blk;
    //     blk.blk_ptr = &blk_pool[blk_start_off]; 
    //     global_blk_vec.push_back(blk);
    //     cout << "blk_start_off: " << blk_start_off << " total size: " << (blk_num*PM_BLK_SIZE) << " malloc_cnt: " << malloc_cnt << endl;

    //     // assert (blk_start_off <= blk_num*PM_BLK_SIZE);
    //     if (memkind_detect_kind(&blk_pool[blk_start_off]) == MEMKIND_DEFAULT) {
    //         cout << "is dram" << endl;
    //     }
    //     // cout << "generate test. blk[0]: " << blk.blk_ptr[0] << endl;
    //     // initializeBitmap(&blk.myBitmap);
    //     blk_start_off += PM_BLK_SIZE;
    // }
    // while(1);
}

void generate_blk_batch_for_parallel(int blk_num, int thread_num, memkind_t kind)
{
    // cout << "generate blk batch!" << endl;
    // blk_num = 10000;
    int32_t *blk_pool = (int32_t *)memkind_malloc(kind, PM_BLK_SIZE * sizeof(int32_t) * blk_num);
    memset(blk_pool, -1, PM_BLK_SIZE * sizeof(int32_t) * blk_num);
    // global_array_vec.push_back(blk_pool);
    // cout << "blk_pool[0]: " << blk_pool[1000] << endl;
    int32_t blk_start_off = 0;
    for (int i = 0; i < blk_num; i++) {
        pm_blk blk;
        blk.blk_ptr = &blk_pool[blk_start_off]; 
        global_parallel_blk_vec[thread_num].push_back(blk);
        // cout << "generate test. blk[0]: " << blk.blk_ptr[0] << endl;
        // initializeBitmap(&blk.myBitmap);
        blk_start_off += PM_BLK_SIZE;
    }
}

PMA* pma_creat (memkind_t kind)
{
    PMA* pma = (PMA*)malloc(sizeof(PMA));
    pma->ele_num = 0;
    pma->arr_cap = 4;
    pma->seg_len = 2;
    pma->seg_num = pma->arr_cap / pma->seg_len;
    pma->tree_h = floor_lg (pma->seg_num) + 1;
    pma->delta_t = (t_0 - t_h) / pma->tree_h;  
    pma->delta_p = (p_h - p_0) / pma->tree_h;
    // cout << "tree_h: " << pma->tree_h << " seg_num: " << pma->seg_num << endl;
    // pma->array = (element_type *)memkind_malloc(kind, sizeof(element_type) * pma->arr_cap);  // on NVM
    // pma->row_offset = (vertexOff *)malloc(sizeof(vertexOff) * pma->vertex_cap); // on DRAM
    return (pma);
}

int32_t blk_search(pm_blk blk, int32_t des)
{ 
    int32_t index = -1;
    int32_t from = 0;
    int32_t to = PM_BLK_SIZE - 1;

    while (from <= to) {
        int32_t mid = from + (to - from)/2;
        int32_t i = mid;
        // Start scanning left until we find a non-empty slot or we reach past the
        while (i>=from && blk.blk_ptr[i]<0) {
            i --;
        }

        if (i < from) {
            from = mid + 1;
        } else {
            if (blk.blk_ptr[i] == des) {
                index = i;
                break;
            } else if (blk.blk_ptr[i] < des) {
                from = mid + 1;
            } else {
                to = i - 1;
            }
        }
    }
    index = to;
    while (index >= 0 && blk.blk_ptr[index] < 0){
        index --;
    }
    
    return index;
}

int32_t get_blk_min(pm_blk blk)
{
    // TBD: prefetch
    for (int i=0; i<PM_BLK_SIZE; i++) {
        if (blk.blk_ptr[i] >= 0) {
            return blk.blk_ptr[i];
        }
    }
    return -1;
}

int32_t get_blk_max(pm_blk blk)
{
    // TBD: prefetch
    for (int i=PM_BLK_SIZE-1; i>=0; i--) {
        if (blk.blk_ptr[i] >= 0) {
            return blk.blk_ptr[i];
        }
    }
    return -1;
}


int32_t pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from_blk, int32_t to_blk, int32_t &index) 
{
    
    if (src->blk_list.size() == 1)  
    {
        // cout << "test---" << endl;
        index = blk_search(src->blk_list[0], des);
        return 0; 
    }

    int32_t init_to_blk = to_blk;
    index = 0;

    while (from_blk <= to_blk)
    {
        int32_t mid_blk = from_blk + (to_blk - from_blk) / 2;
        pm_blk blk = src->blk_list[mid_blk];
        int32_t cur_blk_max = get_blk_max(blk);
        int32_t cur_blk_min = get_blk_min(blk);
        // cout << "cur_blk_max: " << cur_blk_max << " cur_blk_min: " << cur_blk_min << endl;

        if (des > cur_blk_min && des < cur_blk_max) 
        {
            // cout << "test---" << endl;
            index = blk_search(blk, des);
            return mid_blk;  
        } else if (des < cur_blk_min) {
            if (mid_blk == 0) {  
                index = -1;
                return 0;
            }
            to_blk = mid_blk - 1;
        } else { // 大于blk.max_nei
            if (mid_blk != init_to_blk) { 
                // int32_t next_min_nei = src.blk_list[mid_blk+1].min_nei;
                int32_t next_min_nei = get_blk_min(src->blk_list[mid_blk+1]);
                if (des < next_min_nei) { 
                     index = blk_search(blk, des); 
                     return mid_blk;
                }
            }
            if (from_blk == to_blk) {
                index = blk_search(blk, des); 
                return mid_blk;
            }
            from_blk = mid_blk + 1;
        }
    }
}

int32_t binary_pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to) {
    int32_t index;
    int32_t blk_idx, bit_idx;
    while (from <= to) {
        int32_t mid = from + (to - from) / 2;
        int32_t i = mid;
        /* Start scanning left until we find a non-empty slot or we reach past the
        * beginning of the subarray. */
        blk_idx = i / PM_BLK_SIZE;
        bit_idx = i % PM_BLK_SIZE;
        while (i >= from && (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)){
            i--;
            blk_idx = i / PM_BLK_SIZE;
            bit_idx = i % PM_BLK_SIZE;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx--;
            // }
        }

        if (i < from) {  /* Everything between from and mid (inclusive) is empty. */
            from = mid + 1;
        } else {
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] == des) {
                index = i;
                // cout << "error: insert the same neighbor. src: " << src << " des: " << des << endl;
                return index;
            }
            else if (src->blk_list[blk_idx].blk_ptr[bit_idx] < des) {
                from = mid + 1;
            } else {  /* pma->array [i].key > key */
                to = i - 1;  
            }
        }
    }
    /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
    index = to;
    blk_idx = index / PM_BLK_SIZE;
    bit_idx = index % PM_BLK_SIZE;
    while (index >= 0 && (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)) {
        index--;   
        blk_idx = index / PM_BLK_SIZE;
        bit_idx = index % PM_BLK_SIZE;
        // bit_idx--;
        // if (bit_idx < 0) {
        //     bit_idx = PM_BLK_SIZE - 1;
        //     blk_idx--;
        // }
    }
    return index;
}

int32_t binary_pma_find_mix(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to) {
    int32_t index;
    int32_t blk_idx, bit_idx;
    while (from <= to) {
        int32_t mid = from + (to - from) / 2;
        int32_t i = mid;
        /* Start scanning left until we find a non-empty slot or we reach past the
        * beginning of the subarray. */
        blk_idx = i / PM_BLK_SIZE;
        bit_idx = i % PM_BLK_SIZE;
        while (i >= from && (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)){
            i--;
            blk_idx = i / PM_BLK_SIZE;
            bit_idx = i % PM_BLK_SIZE;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx--;
            // }
        }

        if (i < from) {  /* Everything between from and mid (inclusive) is empty. */
            from = mid + 1;
        } else {
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] == des) {
                index = i;
                
                src->blk_list[blk_idx].blk_ptr[bit_idx] = -1;
                src->degree --;
                pma->ele_num--;
                // cout << "error: insert the same neighbor. src: " << src << " des: " << des << endl;
                return index;
            }
            else if (src->blk_list[blk_idx].blk_ptr[bit_idx] < des) {
                from = mid + 1;
            } else {  /* pma->array [i].key > key */
                to = i - 1;  
            }
        }
    }
    /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
    index = to;
    blk_idx = index / PM_BLK_SIZE;
    bit_idx = index % PM_BLK_SIZE;
    while (index >= 0 && (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)) {
        index--;  
        blk_idx = index / PM_BLK_SIZE;
        bit_idx = index % PM_BLK_SIZE;
        // bit_idx--;
        // if (bit_idx < 0) {
        //     bit_idx = PM_BLK_SIZE - 1;
        //     blk_idx--;
        // }
    }
    return index;
}

int32_t pma_insert(PMA* pma, Vertex *src, int32_t des, int32_t index, memkind_t kind)
{
    /*
        insert after index 
    */
    int32_t j = index + 1;
    int32_t blk_idx, bit_idx;

    // cout << "j: " << j << " pma->arr_cap: " <<  pma->arr_cap << endl;
    // blk_idx = j / PM_BLK_SIZE;  
    // bit_idx = j % PM_BLK_SIZE;
    while (j < pma->arr_cap) { 
        blk_idx = j / PM_BLK_SIZE;  
        bit_idx = j % PM_BLK_SIZE;
        // pm_blk blk = src.blk_list[blk_idx];
        if (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)
            break;
        // if (!testBit(&blk.myBitmap, bit_idx))
        //     break;
        j++;
        // bit_idx++;
        // if (bit_idx == PM_BLK_SIZE) {
        //     blk_idx ++;
        //     bit_idx = 0;
        // }
    }
    
    if (j < pma->arr_cap) {
        // cout << "test --0 des: " << des << endl;
        // blk_idx = j / PM_BLK_SIZE;  
        // bit_idx = j % PM_BLK_SIZE;
        // int32_t pre_blk_idx = (j-1) / PM_BLK_SIZE;  
        // int32_t pre_bit_idx = (j-1) % PM_BLK_SIZE;

        while (j > index + 1) /* Push elements to make space for the new element. */
        {
            blk_idx = j / PM_BLK_SIZE;  
            bit_idx = j % PM_BLK_SIZE;
            int32_t pre_blk_idx = (j-1) / PM_BLK_SIZE;  
            int32_t pre_bit_idx = (j-1) % PM_BLK_SIZE;

            src->blk_list[blk_idx].blk_ptr[bit_idx] = src->blk_list[pre_blk_idx].blk_ptr[pre_bit_idx];

            j--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx --;
            // }

            // pre_bit_idx--;
            // if (pre_bit_idx < 0) {
            //     pre_bit_idx = PM_BLK_SIZE - 1;
            //     pre_blk_idx --;
            // }
        }

        blk_idx = (index + 1) / PM_BLK_SIZE;  
        bit_idx = (index + 1) % PM_BLK_SIZE;
        src->blk_list[blk_idx].blk_ptr[bit_idx] = des;
        // cout << "test--- bit_idx: " << bit_idx << endl;
        index = index + 1;  
    } else { /* No empty space to the right of i. Try left. */
        // cout << "test --1 des: " << des << " index: " << index << endl;
        j = index - 1;
        // blk_idx = j / PM_BLK_SIZE;  
        // bit_idx = j % PM_BLK_SIZE;
        while (j > 0) {
            blk_idx = j / PM_BLK_SIZE;  
            bit_idx = j % PM_BLK_SIZE;
            pm_blk blk = src->blk_list[blk_idx];
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)
                break;
            // if (!testBit(&blk.myBitmap, bit_idx))
            //     break;
            j--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx --;
            // }
        }
        if (j >= 0) {
            // blk_idx = j / PM_BLK_SIZE;  
            // bit_idx = j % PM_BLK_SIZE;
            // int32_t be_blk_idx = (j+1) / PM_BLK_SIZE;  
            // int32_t be_bit_idx = (j+1) % PM_BLK_SIZE;
            // pm_blk blk = src.blk_list[blk_idx];
            // setBit(&src.blk_list[blk_idx].myBitmap, bit_idx);
            while (j < index) {
                blk_idx = j / PM_BLK_SIZE;  
                bit_idx = j % PM_BLK_SIZE;
                int32_t be_blk_idx = (j+1) / PM_BLK_SIZE;  
                int32_t be_bit_idx = (j+1) % PM_BLK_SIZE;
                src->blk_list[blk_idx].blk_ptr[bit_idx] = src->blk_list[be_blk_idx].blk_ptr[be_bit_idx];
                j++;

                // bit_idx++;
                // if (bit_idx == PM_BLK_SIZE) {
                //     bit_idx = 0;
                //     blk_idx++;
                // }
                // be_bit_idx++;
                // if (be_bit_idx == PM_BLK_SIZE) {
                //     be_bit_idx = 0;
                //     be_blk_idx++;
                // }
            }
            blk_idx = (index) / PM_BLK_SIZE;  
            bit_idx = (index) % PM_BLK_SIZE;
            src->blk_list[blk_idx].blk_ptr[bit_idx] = des;
        }
    }
    pma->ele_num++;
    src->degree++;
    // print_graph();

    // cout << "test---" << endl;
    
    rebalance(pma, src, index, kind);

    return index;
}

int32_t pma_insert_mv(PMA* pma, Vertex *src, int32_t des, int32_t index, memkind_t kind)
{
    /*
        insert after index 
        insert the new des in current blk
    */

    int32_t j = index + 1;
    int32_t blk_idx, bit_idx;

    // cout << "j: " << j << " pma->arr_cap: " <<  pma->arr_cap << endl;
    // blk_idx = j / PM_BLK_SIZE;  
    // bit_idx = j % PM_BLK_SIZE;
    while (j < pma->arr_cap) { 
        blk_idx = j / PM_BLK_SIZE;  
        bit_idx = j % PM_BLK_SIZE;
        // pm_blk blk = src.blk_list[blk_idx];
        if (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)
            break;
        // if (!testBit(&blk.myBitmap, bit_idx))
        //     break;
        j++; 
        // bit_idx++;
        // if (bit_idx == PM_BLK_SIZE) {
        //     blk_idx ++;
        //     bit_idx = 0;
        // }
    }
    
    if (j < pma->arr_cap) {
        // cout << "insert right..." << endl;
        int32_t max_blk = j / PM_BLK_SIZE;
        int32_t min_blk = (index + 1) / PM_BLK_SIZE;
        pm_blk* next_blk = NULL;
        pm_blk* cur_blk = NULL;
        // cout << "min_blk: " << min_blk << " max_blk: " << max_blk << endl;
        for (int32_t tmp_blk = min_blk; tmp_blk <= max_blk; tmp_blk++) {
            if (src->blk_list[tmp_blk].next_blk == NULL) {
                // cout << " test --- 0 global_task_id: " << global_task_id << endl;
                // pm_blk blk = generate_blk(kind);
                pm_blk* blk = new pm_blk(generate_blk(kind));
                // blk.next_blk = NULL;
                // if (blk.next_blk == NULL)
                //     cout << "test --- 0.5 blk.nextblk is NULL" << endl;
                blk->task_id = global_task_id;
                memcpy(blk->blk_ptr, src->blk_list[tmp_blk].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
                src->blk_list[tmp_blk].next_blk = blk;
                // src->blk_list[tmp_blk].next_blk->task_id = 2;
                // if (src->blk_list[tmp_blk].next_blk->next_blk == NULL)
                //     cout << "test --- 1 last blk task id: " << src->blk_list[tmp_blk].next_blk->task_id << " real: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
            } else {
                cur_blk = &(src->blk_list[tmp_blk]);
                while (cur_blk->next_blk != NULL) {
                    cur_blk = cur_blk->next_blk;
                }
                if (cur_blk->task_id < global_task_id) {
                    // pm_blk blk = generate_blk(kind);
                    pm_blk* blk = new pm_blk(generate_blk(kind));
                    blk->task_id = global_task_id;
                    memcpy(blk->blk_ptr, src->blk_list[tmp_blk].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
                    cur_blk->next_blk = blk;
                }
            }
        }

        while (j > index + 1) /* Push elements to make space for the new element. */
        {
            blk_idx = j / PM_BLK_SIZE;  
            bit_idx = j % PM_BLK_SIZE;
            int32_t pre_blk_idx = (j-1) / PM_BLK_SIZE;  
            int32_t pre_bit_idx = (j-1) % PM_BLK_SIZE;

            // if (src->blk_list[blk_idx].next_blk == NULL) {
            //     if (src->blk_list[blk_idx].task_id < global_task_id) {
            //         pm_blk blk = generate_blk(kind);  
            //         blk.task_id = global_task_id;
            //         memcpy(blk.blk_ptr, src->blk_list[blk_idx].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
            //         src->blk_list[blk_idx].next_blk = &blk;
            //     }
            // } else {
            //     pm_blk* next_blk = src->blk_list[blk_idx].next_blk; // 
            //     while (next_blk->next_blk != NULL) {
            //         next_blk = next_blk->next_blk;
            //     }
            //     if (next_blk->task_id < global_task_id) {
            //         pm_blk blk = generate_blk(kind);  // 
            //         blk.task_id = global_task_id;
            //         next_blk->next_blk = &blk;
            //         memcpy(blk.blk_ptr, src->blk_list[blk_idx].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
            //     }
            // }

            // if (src->blk_list[blk_idx].task_id < global_task_id && src->blk_list[blk_idx].next_blk == NULL) 
            // {
            //     pm_blk blk = generate_blk(kind);  // 
            //     blk.task_id = global_task_id;
            //     memcpy(blk.blk_ptr, src->blk_list[blk_idx].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
            //     src->blk_list[blk_idx].next_blk = &blk;
            // } else if (src->blk_list[blk_idx].next_blk != NULL) {
                
            // }


            src->blk_list[blk_idx].blk_ptr[bit_idx] = src->blk_list[pre_blk_idx].blk_ptr[pre_bit_idx];

            j--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx --;
            // }

            // pre_bit_idx--;
            // if (pre_bit_idx < 0) {
            //     pre_bit_idx = PM_BLK_SIZE - 1;
            //     pre_blk_idx --;
            // }
        }

        blk_idx = (index + 1) / PM_BLK_SIZE;  // 
        bit_idx = (index + 1) % PM_BLK_SIZE;
        src->blk_list[blk_idx].blk_ptr[bit_idx] = des;
        // cout << "test--- bit_idx: " << bit_idx << endl;
        index = index + 1;  
    } else { /* No empty space to the right of i. Try left. */
        // cout << "test --1 des: " << des << " index: " << index << endl;
        // cout << "insert left..." << endl;
        j = index - 1;
        // blk_idx = j / PM_BLK_SIZE;  
        // bit_idx = j % PM_BLK_SIZE;
        while (j > 0) {
            blk_idx = j / PM_BLK_SIZE;  
            bit_idx = j % PM_BLK_SIZE;
            pm_blk blk = src->blk_list[blk_idx];
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)
                break;
            // if (!testBit(&blk.myBitmap, bit_idx))
            //     break;
            j--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx --;
            // }
        }
        if (j >= 0) {
            int32_t max_blk = index / PM_BLK_SIZE;
            int32_t min_blk = j / PM_BLK_SIZE;
            pm_blk* next_blk = NULL;
            pm_blk* cur_blk = NULL;
            for (int32_t tmp_blk = min_blk; tmp_blk <= max_blk; tmp_blk++) {
                if (src->blk_list[tmp_blk].next_blk == NULL) {
                    // pm_blk blk = generate_blk(kind);
                    pm_blk* blk = new pm_blk(generate_blk(kind));
                    blk->task_id = global_task_id;
                    memcpy(blk->blk_ptr, src->blk_list[tmp_blk].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
                    src->blk_list[tmp_blk].next_blk = blk;
                } else {
                    cur_blk = &(src->blk_list[tmp_blk]);
                    while (cur_blk->next_blk != NULL) {
                        cur_blk = cur_blk->next_blk;
                    }
                    if (cur_blk->task_id < global_task_id) {
                        // pm_blk blk = generate_blk(kind);
                        pm_blk* blk = new pm_blk(generate_blk(kind));
                        blk->task_id = global_task_id;
                        memcpy(blk->blk_ptr, src->blk_list[tmp_blk].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
                        cur_blk->next_blk = blk;
                    }
                }
            }
            // blk_idx = j / PM_BLK_SIZE;  
            // bit_idx = j % PM_BLK_SIZE;
            // int32_t be_blk_idx = (j+1) / PM_BLK_SIZE;  
            // int32_t be_bit_idx = (j+1) % PM_BLK_SIZE;
            // pm_blk blk = src.blk_list[blk_idx];
            // setBit(&src.blk_list[blk_idx].myBitmap, bit_idx);
            while (j < index) {
                blk_idx = j / PM_BLK_SIZE;  
                bit_idx = j % PM_BLK_SIZE;
                int32_t be_blk_idx = (j+1) / PM_BLK_SIZE;  
                int32_t be_bit_idx = (j+1) % PM_BLK_SIZE;
                src->blk_list[blk_idx].blk_ptr[bit_idx] = src->blk_list[be_blk_idx].blk_ptr[be_bit_idx];
                j++;

                // bit_idx++;
                // if (bit_idx == PM_BLK_SIZE) {
                //     bit_idx = 0;
                //     blk_idx++;
                // }
                // be_bit_idx++;
                // if (be_bit_idx == PM_BLK_SIZE) {
                //     be_bit_idx = 0;
                //     be_blk_idx++;
                // }
            }
            blk_idx = (index) / PM_BLK_SIZE;  
            bit_idx = (index) % PM_BLK_SIZE;
            src->blk_list[blk_idx].blk_ptr[bit_idx] = des;
        }
    }
    pma->ele_num++;
    src->degree++;
    // print_graph();


    // cout << "test---" << endl;
    // cout << "v1 next blk task id: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;

    
    rebalance_mv(pma, src, index, kind);

    // cout << "v1 next blk task id 0.1: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
    // cout << "v1 next blk task id 0.2: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
    // cout << "v1 next blk task id 0.3: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
    return index;
}


void rebalance(PMA* pma, Vertex *src, int32_t index, memkind_t kind) 
{
    int32_t window_start, window_end;
    int32_t height = 0;
    uint32_t occupancy = 1;
    int32_t left_index = index - 1;
    int32_t right_index = index + 1;
    double density, t_height, p_height;
    int32_t blk_idx, bit_idx;
    
    do {
        uint32_t window_size = pma->seg_len * (1 << height);
        uint32_t window = index / window_size;
        window_start = window * window_size;
        window_end = window_start + window_size;
        // blk_idx = left_index / PM_BLK_SIZE;
        // bit_idx = left_index % PM_BLK_SIZE;
        while (left_index >= window_start) {
            blk_idx = left_index / PM_BLK_SIZE;
            bit_idx = left_index % PM_BLK_SIZE;
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] >= 0)
                occupancy++;
            // if (testBit(&src.blk_list[blk_idx].myBitmap, bit_idx))
            //     occupancy++;
            left_index--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx--;
            // }
        }
        while (right_index < window_end) {
            blk_idx = right_index / PM_BLK_SIZE;
            bit_idx = right_index % PM_BLK_SIZE;
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] >= 0)
                occupancy++;
            // if (testBit(&src.blk_list[blk_idx].myBitmap, bit_idx))
            //     occupancy++;
            right_index++;
            // bit_idx++;
            // if (bit_idx == PM_BLK_SIZE) {
            //     bit_idx = 0;
            //     blk_idx++;
            // }
        }
        density = (double)occupancy / (double)window_size;
        t_height = t_0 - (height * pma->delta_t);
        p_height = p_0 + (height * pma->delta_p);
        height++;
    } while (( 
            density >= t_height) &&
           height < pma->tree_h);
    // while ((density < p_height ||  //  density < p_height
    //         density >= t_height) &&
    //        height < pma->tree_h);

    // cout << "height: " << height << " pma tree: " << pma->tree_h << endl;
    // cout << "occupancy: " << occupancy << " window_start: " << window_start << " window_end: " << window_end << endl;
    // print_pma(src);
    // print_vertex(src);
      /* Found a window within threshold. */
    if (height == 1)
        return;

    // cout << "test----" << endl;
    // int32_t tmpOccupacy = 0;
    // for (int32_t i = window_start; i < window_end; i++) {
    //     int32_t tmpBlkIdx = i / PM_BLK_SIZE;
    //     int32_t tmpBitIdx = i % PM_BLK_SIZE;
    //     if (src->blk_list[tmpBlkIdx].blk_ptr[tmpBitIdx] >= 0) {
    //         tmpOccupacy++;
    //     }
    // }
    // if (tmpOccupacy != occupancy)
    //     cout << "tmpOccupacy: " << tmpOccupacy << " occupancy: " << occupancy << endl;

    
    if (density >= p_height && density < t_height) {
        // pack (pma, window_start, window_end, occupancy);
        // cout << "test --- rebalance spread. occupancy: " << occupancy << endl;
        // print_pma(src);
        // print_vertex(src);
        spread (pma, src, window_start, window_end, occupancy, kind);
        // print_vertex(src);
    } else {
        // cout << "test --- rebalance resize. occupancy: " << occupancy << endl;
        resize (pma, src, kind);
        // print_pma(src);
        // print_vertex(src);
  }
}

void rebalance_mv(PMA* pma, Vertex *src, int32_t index, memkind_t kind) 
{
    int32_t window_start, window_end;
    int32_t height = 0;
    uint32_t occupancy = 1;
    int32_t left_index = index - 1;
    int32_t right_index = index + 1;
    double density, t_height, p_height;
    int32_t blk_idx, bit_idx;
    
    do {
        uint32_t window_size = pma->seg_len * (1 << height);
        uint32_t window = index / window_size;
        window_start = window * window_size;
        window_end = window_start + window_size;
        // blk_idx = left_index / PM_BLK_SIZE;
        // bit_idx = left_index % PM_BLK_SIZE;
        while (left_index >= window_start) {
            blk_idx = left_index / PM_BLK_SIZE;
            bit_idx = left_index % PM_BLK_SIZE;
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] >= 0)
                occupancy++;
            // if (testBit(&src.blk_list[blk_idx].myBitmap, bit_idx))
            //     occupancy++;
            left_index--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx--;
            // }
        }
        while (right_index < window_end) {
            blk_idx = right_index / PM_BLK_SIZE;
            bit_idx = right_index % PM_BLK_SIZE;
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] >= 0)
                occupancy++;
            // if (testBit(&src.blk_list[blk_idx].myBitmap, bit_idx))
            //     occupancy++;
            right_index++;
            // bit_idx++;
            // if (bit_idx == PM_BLK_SIZE) {
            //     bit_idx = 0;
            //     blk_idx++;
            // }
        }
        density = (double)occupancy / (double)window_size;
        t_height = t_0 - (height * pma->delta_t);
        p_height = p_0 + (height * pma->delta_p);
        height++;
    } while (density >= t_height && height < pma->tree_h);  
           

    // cout << "height: " << height << " pma tree: " << pma->tree_h << endl;
    // cout << "occupancy: " << occupancy << " window_start: " << window_start << " window_end: " << window_end << endl;
    // print_pma(src);
    // print_vertex(src);
      /* Found a window within threshold. */
    if (height == 1)
        return;

    // cout << "rebalance test----" << endl;
    // int32_t tmpOccupacy = 0;
    // for (int32_t i = window_start; i < window_end; i++) {
    //     int32_t tmpBlkIdx = i / PM_BLK_SIZE;
    //     int32_t tmpBitIdx = i % PM_BLK_SIZE;
    //     if (src->blk_list[tmpBlkIdx].blk_ptr[tmpBitIdx] >= 0) {
    //         tmpOccupacy++;
    //     }
    // }
    // if (tmpOccupacy != occupancy)
    //     cout << "tmpOccupacy: " << tmpOccupacy << " occupancy: " << occupancy << endl;

    
    if (density >= p_height && density < t_height) {
        int32_t max_blk = (window_end-1) / PM_BLK_SIZE;
        int32_t min_blk = window_start / PM_BLK_SIZE;
        pm_blk* next_blk = NULL;
        pm_blk* cur_blk = NULL;
        for (int32_t tmp_blk = min_blk; tmp_blk <= max_blk; tmp_blk++) {
            if (src->blk_list[tmp_blk].next_blk == NULL) {
                // pm_blk blk = generate_blk(kind);
                pm_blk* blk = new pm_blk(generate_blk(kind));
                blk->task_id = global_task_id;
                memcpy(blk->blk_ptr, src->blk_list[tmp_blk].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
                src->blk_list[tmp_blk].next_blk = blk;
            } else {
                cur_blk = &(src->blk_list[tmp_blk]);
                while (cur_blk->next_blk != NULL) {
                    cur_blk = cur_blk->next_blk;
                }
                if (cur_blk->task_id < global_task_id) {
                    // pm_blk blk = generate_blk(kind);
                    pm_blk* blk = new pm_blk(generate_blk(kind));
                    blk->task_id = global_task_id;
                    memcpy(blk->blk_ptr, src->blk_list[tmp_blk].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
                    cur_blk->next_blk = blk;
                }
            }
        }
        // pack (pma, window_start, window_end, occupancy);
        // cout << "test --- rebalance spread. occupancy: " << occupancy << endl;
        // print_pma(src);
        // print_vertex(src);
        spread (pma, src, window_start, window_end, occupancy, kind);
        // print_vertex(src);
    } else {
        // cout << "test --- rebalance resize. occupancy: " << occupancy << endl;
        int32_t max_blk = src->blk_list.size();
        int32_t min_blk = 0;
        pm_blk* next_blk = NULL;
        pm_blk* cur_blk = NULL;
        for (int32_t tmp_blk = min_blk; tmp_blk < max_blk; tmp_blk++) {
            if (src->blk_list[tmp_blk].next_blk == NULL) {
                // pm_blk blk = generate_blk(kind);
                pm_blk* blk = new pm_blk(generate_blk(kind));
                blk->task_id = global_task_id;
                memcpy(blk->blk_ptr, src->blk_list[tmp_blk].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
                src->blk_list[tmp_blk].next_blk = blk;
            } else {
                cur_blk = &(src->blk_list[tmp_blk]);
                while (cur_blk->next_blk != NULL) {
                    cur_blk = cur_blk->next_blk;
                }
                if (cur_blk->task_id < global_task_id) {
                    // pm_blk blk = generate_blk(kind);
                    pm_blk* blk = new pm_blk(generate_blk(kind));
                    blk->task_id = global_task_id;
                    memcpy(blk->blk_ptr, src->blk_list[tmp_blk].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
                    cur_blk->next_blk = blk;
                }
            }
        }
        resize (pma, src, kind);
        for (int32_t tmp_blk = max_blk; tmp_blk < src->blk_list.size(); tmp_blk++) {
            src->blk_list[tmp_blk].task_id = global_task_id;
        }
        // print_pma(src);
        // print_vertex(src);
  }
}

void spread (PMA* pma, Vertex *src, int32_t from, int32_t to, int32_t occupancy, memkind_t kind)
{
 
    assert (from < to);
    int32_t capacity = to - from; 
    int32_t frequency = (capacity << 8) / occupancy;  /* 8-bit fixed point arithmetic. */
    int32_t read_index = to - 1;  
    int32_t write_index = (to << 8) - frequency;

    int32_t start_blk = from / PM_BLK_SIZE;
    int32_t end_blk = (to - 1) / PM_BLK_SIZE;

    int32_t span_blks = end_blk - start_blk + 1;

    
    int32_t *tmp_win_arr = (int32_t*)malloc(capacity * sizeof(int32_t));
    // int32_t *tmp_win_arr = (int32_t*)memkind_malloc(kind, capacity * sizeof(int32_t));

    // memset(tmp_win_arr, -1, capacity * sizeof(int32_t));

    if (tmp_win_arr != NULL) {
    //  -1
        for (size_t i = 0; i < capacity; ++i) {
            tmp_win_arr[i] = -1;
        }
    }

    // test
    // cout << "test---- span_blks: " << span_blks << " from: " << from << " to: " << to
    //     << " start_blk: " << start_blk << " end_blk: " << end_blk << endl;
    // int32_t tmpOccupacy = 0;
    // for (int32_t i = from; i < to; i++) {
    //     int32_t tmpBlkIdx = i / PM_BLK_SIZE;
    //     int32_t tmpBitIdx = i % PM_BLK_SIZE;
    //     if (src->blk_list[tmpBlkIdx].blk_ptr[tmpBitIdx] >= 0) {
    //         tmpOccupacy++;
    //     }
    // }
    // cout << "tmpOccupacy: " << tmpOccupacy << " occupancy: " << occupancy << endl;

    // cout << "span_blks: " << span_blks << endl;
    // cout << "from: " << from << " to: " << to << endl;
    // print_vertex(src);
    if (span_blks == 1) {
        int32_t base_off = start_blk * PM_BLK_SIZE;
        read_index = read_index - base_off;
        write_index = ((to-base_off) << 8) - frequency;  
        from = from - base_off;
        // cout << "from: " << from << " to: " << to << " write_index: " << (write_index >> 8) << endl;
        // Bitmap32 readBitmap = src.blk_list[start_blk].myBitmap; 
        // print_bitmap(&readBitmap);
        //  bitmap
        // for (int i = from; i < (to-base_off); i++) {
        //     clearBit(&src.blk_list[start_blk].myBitmap, i);
        // }
        int32_t real_write_idx = (write_index >> 8) - from;
        occupancy -= 1; 
        while (read_index >= from) {
            if (src->blk_list[start_blk].blk_ptr[read_index] >= 0) {
            // if (testBit(&readBitmap, read_index)) {
                // src.blk_list[0].blk_ptr[write_index >> 8] = src.blk_list[0].blk_ptr[read_index];
                assert(real_write_idx >= 0);
                tmp_win_arr[real_write_idx] = src->blk_list[start_blk].blk_ptr[read_index];
                // clearBit(&src.blk_list[start_blk].myBitmap, read_index);
                // cout << "write_index " << (write_index >> 8) - from << " read_index: " << read_index << " des: " << src.blk_list[start_blk].blk_ptr[read_index] << endl;
                // setBit(&src.blk_list[start_blk].myBitmap, (write_index >> 8));
                write_index -= frequency;
                read_index--;
                occupancy --;
                real_write_idx = (write_index >> 8) - from;
                if (real_write_idx < occupancy) {
                    real_write_idx = occupancy;
                }

            } else {
                read_index --;
            }
        }
        // copy tmp_win_arr to pma
        memcpy(&src->blk_list[start_blk].blk_ptr[from], tmp_win_arr, sizeof(int32_t) * capacity); 
        // free(tmp_win_arr);
        // if (from < 0 || from >= PM_BLK_SIZE || (from + capacity - 1) >= 32) {
        //     cout << "error 1" << endl;
        //     exit(1);
        // }
    } else {  
        // write_index = (capacity << 8) - frequency; 
        // vector<Bitmap32> resultBitMapVec;
        // for (int32_t j = start_blk; j <= end_blk; j++) {
        //     Bitmap32 bitMap;
        //     initializeBitmap(&bitMap);
        //     resultBitMapVec.push_back(bitMap);
        //     
        // }
        // int32_t tmpResBitmapIdx = span_blks - 1; 
        // write_index = (occupancy << 8) - frequency;
        occupancy -= 1; 
        int32_t real_write_idx = (write_index >> 8) - from;
        for (int i = end_blk; i >= start_blk; i--) {
            int32_t loc_win_from, loc_win_to; 
            if (i > start_blk) {
                loc_win_from = 0;
            } else {
                loc_win_from = from - i * PM_BLK_SIZE;
            }
            if (i < end_blk) {
                loc_win_to = PM_BLK_SIZE - 1;
            } else {
                loc_win_to = to - i * PM_BLK_SIZE - 1;
            }
            // cout << "neighbor: ";
            while (loc_win_to >= loc_win_from) {
                if (src->blk_list[i].blk_ptr[loc_win_to] >= 0) {
                    // cout << "tmp_cnt: " << tmp_cnt << endl;
                    assert(occupancy >= 0);
                // if (testBit(&src.blk_list[i].myBitmap, loc_win_to)) {
                    // tmp_win_arr[(write_index >> 8) - from] = src->blk_list[i].blk_ptr[loc_win_to];
                    tmp_win_arr[real_write_idx] = src->blk_list[i].blk_ptr[loc_win_to];
                    // cout << src->blk_list[i].blk_ptr[loc_win_to] << ", ";
                    // clearBit(&src.blk_list[i].myBitmap, loc_win_to);
                    // int32_t blk_offset = (write_index >> 8) % PM_BLK_SIZE;
                    // setBit(&resultBitMapVec[tmpResBitmapIdx], blk_offset); 
                    write_index -= frequency;
                    real_write_idx = (write_index >> 8) - from;
                    occupancy --;

                    if (real_write_idx < occupancy) {
                        real_write_idx = occupancy;
                    }
                    
                }
                loc_win_to--;
            }
            // assert((write_index >> 8) >= from);
        }
        // cout << endl;

        // cout << "tmp_win_arr: ";
        // for (int i = 0; i < capacity; i++) {
        //     cout << tmp_win_arr[i] << ", ";
        // }
        // cout << endl;

        // tmpResBitmapIdx = 0;
        int32_t cpy_src = 0;
        int32_t cpy_des = from % PM_BLK_SIZE;
        int32_t cpy_size = PM_BLK_SIZE - cpy_des;   
        // cout << "start blk: " << start_blk << " end blk: " << end_blk << endl; 
        for (int32_t i = start_blk; i <= end_blk; i++) {
            // cout << "cpy_des: " << cpy_des << " cpy_src: " << cpy_src << " cpy_size: " << cpy_size << endl;
            // andBit(&src.blk_list[i].myBitmap, &resultBitMapVec[tmpResBitmapIdx]);
            // cout << "bitmap des: " << src.blk_list[i].myBitmap.word << " bitmap src: " << resultBitMapVec[tmpResBitmapIdx].word << endl;
            // tmpResBitmapIdx++;
            memcpy(&src->blk_list[i].blk_ptr[cpy_des], &tmp_win_arr[cpy_src], cpy_size * sizeof(int32_t));
            
            assert((cpy_des+cpy_size-1)<PM_BLK_SIZE);
            
            cpy_des = 0;
            cpy_src += cpy_size;

            if (i == (end_blk-1)) { 
                cpy_size = ((to-1) % PM_BLK_SIZE) + 1;  
            } else {
                cpy_size = PM_BLK_SIZE;
            }
        }
        // free(tmp_win_arr);
        // print_vertex(src);
    }
    free(tmp_win_arr);
}

void resize (PMA* pma, Vertex *src, memkind_t kind)
{
    int32_t read_index = pma->arr_cap - 1;  
    int32_t readBitmapIdx = src->blk_list.size() - 1;
    // print_vertex(src);
    compute_capacity(pma);
    pma->tree_h = floor_lg (pma->seg_num) + 1;
    pma->delta_t = (t_0 - t_h) / pma->tree_h;
    pma->delta_p = (p_h - p_0) / pma->tree_h;


    int32_t blk_num = ceil((double)pma->arr_cap / PM_BLK_SIZE);  
    int32_t ex_blk_num = blk_num - src->blk_list.size();


    // cout << "ex_blk_num: " << ex_blk_num << endl;
    if (global_blk_vec.size() < ex_blk_num) {
        // cout << "ex_blk_num: " << ex_blk_num << endl;
        #pragma omp critical
        {
            generate_blk_batch(ex_blk_num * 2, kind);
        }
    }

    // cout << "blk_num: " << blk_num << " ex_blk_num: " << ex_blk_num << endl;
    // print_vertex(src);
    // print_pma(src);
    // cout << "test -- 0" << endl;

    if (ex_blk_num > 0)
    {
        // cout << "test ---0 id: " << src->id << " degree: " << src->degree << " blk_num: " << blk_num << endl;
        // print_pma(src);
        // src->blk_list = (pm_blk *)realloc(src->blk_list, sizeof(pm_blk) * blk_num);
        // cout << "test ---1" << endl;
        // int32_t src_blk_idx = src->blk_list.size();
        // cout << "global_blk_vec size: " << global_blk_vec.size() << " degree: " << src.degree << endl;
        // cout << "global_blk_vec size: " << global_blk_vec.size() << endl;
        #pragma omp critical
        {
            for (int32_t i=0; i < ex_blk_num; i++) {
                // pm_blk blk = global_blk_vec.back();  
                // initializeBitmap(&blk.myBitmap);
                // cout << "blk[0]: " << blk.blk_ptr[0] << endl;
                
                src->blk_list.push_back(global_blk_vec.back());
                global_blk_vec.pop_back();
                // src->blk_list[src_blk_idx] = blk;
                // src_blk_idx++;
            }
        }

        // src.blk_cnt = blk_num;
    } 

    // cout << "test -- 1" << endl;

    // print_vertex(src);
    // print_pma(src);

    int32_t capacity = pma->arr_cap;
    int32_t frequency = (capacity << 8) / pma->ele_num;
    int32_t write_index = (capacity << 8) - frequency;

    // cout << "read_idx: " << read_index << " write_index: " << (write_index >> 8) << endl;
    int32_t writeBitmapIdx = src->blk_list.size() - 1; 
    
    while ((write_index >> 8) > read_index && read_index >= 0)
    {
        // cout << "read_idx: " << read_index << " write_index: " << (write_index >> 8) 
        //     << " readBitmapIdx: " << readBitmapIdx << " writeBitmapIdx: " << writeBitmapIdx << endl;
        int32_t write_off = (write_index >> 8) % PM_BLK_SIZE;
        int32_t read_off = read_index % PM_BLK_SIZE;
        
        // test
        // cout << "read_off: " << read_off << " write_off: " << write_off << endl;
        // src.blk_list[writeBitmapIdx].blk_ptr[write_off] = 1;
        // cout << "test ---  -1 " << endl;

        if (src->blk_list[readBitmapIdx].blk_ptr[read_off] >= 0) {
        // if (testBit(&src.blk_list[readBitmapIdx].myBitmap, read_off)) {
            // cout << "read_off: " << read_off << " write_off: " << write_off << endl;
            src->blk_list[writeBitmapIdx].blk_ptr[write_off] = src->blk_list[readBitmapIdx].blk_ptr[read_off];
            // cout << "test ---0 " << endl;
            // clearBit(&src.blk_list[readBitmapIdx].myBitmap, read_off);
            src->blk_list[readBitmapIdx].blk_ptr[read_off] = -1;
            // setBit(&src.blk_list[writeBitmapIdx].myBitmap, write_off);
            read_index--;
            write_index -= frequency;

            
            if ((read_index % PM_BLK_SIZE) > read_off) {
                readBitmapIdx --;
                // cout << "readBitmapIdx: " << readBitmapIdx << endl;
            }
            if (((write_index >> 8) % PM_BLK_SIZE) > write_off) {
                writeBitmapIdx--;
                // cout << "writeBitmapIdx: " << writeBitmapIdx << endl;
            }
        } else {
            // cout << "test ---1 " << endl;
            read_index--;
            if ((read_index % PM_BLK_SIZE) > read_off) {
                readBitmapIdx --;
                // cout << "readBitmapIdx: " << readBitmapIdx << endl;
            }
        }
    }
}

int32_t pma_insert_parallel(PMA* pma, Vertex *src, int32_t des, int32_t index, int thread_num, memkind_t kind)
{
    /*
        insert after index 
    */
    int32_t j = index + 1;
    int32_t blk_idx, bit_idx;

    // cout << "j: " << j << " pma->arr_cap: " <<  pma->arr_cap << endl;
    // blk_idx = j / PM_BLK_SIZE;  
    // bit_idx = j % PM_BLK_SIZE;
    while (j < pma->arr_cap) { 
        blk_idx = j / PM_BLK_SIZE;  
        bit_idx = j % PM_BLK_SIZE;
        // pm_blk blk = src.blk_list[blk_idx];
        if (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)
            break;
        // if (!testBit(&blk.myBitmap, bit_idx))
        //     break;
        j++;
        // bit_idx++;
        // if (bit_idx == PM_BLK_SIZE) {
        //     blk_idx ++;
        //     bit_idx = 0;
        // }
    }
    
    if (j < pma->arr_cap) {
        // cout << "test --0 des: " << des << endl;
        // blk_idx = j / PM_BLK_SIZE;  
        // bit_idx = j % PM_BLK_SIZE;
        // int32_t pre_blk_idx = (j-1) / PM_BLK_SIZE;  
        // int32_t pre_bit_idx = (j-1) % PM_BLK_SIZE;

        while (j > index + 1) /* Push elements to make space for the new element. */
        {
            blk_idx = j / PM_BLK_SIZE;  
            bit_idx = j % PM_BLK_SIZE;
            int32_t pre_blk_idx = (j-1) / PM_BLK_SIZE;  
            int32_t pre_bit_idx = (j-1) % PM_BLK_SIZE;

            src->blk_list[blk_idx].blk_ptr[bit_idx] = src->blk_list[pre_blk_idx].blk_ptr[pre_bit_idx];

            j--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx --;
            // }

            // pre_bit_idx--;
            // if (pre_bit_idx < 0) {
            //     pre_bit_idx = PM_BLK_SIZE - 1;
            //     pre_blk_idx --;
            // }
        }

        blk_idx = (index + 1) / PM_BLK_SIZE;  
        bit_idx = (index + 1) % PM_BLK_SIZE;
        src->blk_list[blk_idx].blk_ptr[bit_idx] = des;
        // cout << "test--- bit_idx: " << bit_idx << endl;
        index = index + 1;  
    } else { /* No empty space to the right of i. Try left. */
        // cout << "test --1 des: " << des << " index: " << index << endl;
        j = index - 1;
        // blk_idx = j / PM_BLK_SIZE;  
        // bit_idx = j % PM_BLK_SIZE;
        while (j > 0) {
            blk_idx = j / PM_BLK_SIZE;  
            bit_idx = j % PM_BLK_SIZE;
            pm_blk blk = src->blk_list[blk_idx];
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)
                break;
            // if (!testBit(&blk.myBitmap, bit_idx))
            //     break;
            j--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx --;
            // }
        }
        if (j >= 0) {
            // blk_idx = j / PM_BLK_SIZE; 
            // bit_idx = j % PM_BLK_SIZE;
            // int32_t be_blk_idx = (j+1) / PM_BLK_SIZE;  
            // int32_t be_bit_idx = (j+1) % PM_BLK_SIZE;
            // pm_blk blk = src.blk_list[blk_idx];
            // setBit(&src.blk_list[blk_idx].myBitmap, bit_idx);
            while (j < index) {
                blk_idx = j / PM_BLK_SIZE;  
                bit_idx = j % PM_BLK_SIZE;
                int32_t be_blk_idx = (j+1) / PM_BLK_SIZE;  
                int32_t be_bit_idx = (j+1) % PM_BLK_SIZE;
                src->blk_list[blk_idx].blk_ptr[bit_idx] = src->blk_list[be_blk_idx].blk_ptr[be_bit_idx];
                j++;

                // bit_idx++;
                // if (bit_idx == PM_BLK_SIZE) {
                //     bit_idx = 0;
                //     blk_idx++;
                // }
                // be_bit_idx++;
                // if (be_bit_idx == PM_BLK_SIZE) {
                //     be_bit_idx = 0;
                //     be_blk_idx++;
                // }
            }
            blk_idx = (index) / PM_BLK_SIZE;  
            bit_idx = (index) % PM_BLK_SIZE;
            src->blk_list[blk_idx].blk_ptr[bit_idx] = des;
        }
    }
    pma->ele_num++;
    src->degree++;
    // print_graph();

    // cout << "test---" << endl;
    
    rebalance_parallel(pma, src, index, thread_num, kind);

    return index;
}

void rebalance_parallel(PMA* pma, Vertex *src, int32_t index, int thread_num, memkind_t kind) 
{
    int32_t window_start, window_end;
    int32_t height = 0;
    uint32_t occupancy = 1;
    int32_t left_index = index - 1;
    int32_t right_index = index + 1;
    double density, t_height, p_height;
    int32_t blk_idx, bit_idx;
    
    do {
        uint32_t window_size = pma->seg_len * (1 << height);
        uint32_t window = index / window_size;
        window_start = window * window_size;
        window_end = window_start + window_size;
        // blk_idx = left_index / PM_BLK_SIZE;
        // bit_idx = left_index % PM_BLK_SIZE;
        while (left_index >= window_start) {
            blk_idx = left_index / PM_BLK_SIZE;
            bit_idx = left_index % PM_BLK_SIZE;
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] >= 0)
                occupancy++;
            // if (testBit(&src.blk_list[blk_idx].myBitmap, bit_idx))
            //     occupancy++;
            left_index--;
            // bit_idx--;
            // if (bit_idx < 0) {
            //     bit_idx = PM_BLK_SIZE - 1;
            //     blk_idx--;
            // }
        }
        while (right_index < window_end) {
            blk_idx = right_index / PM_BLK_SIZE;
            bit_idx = right_index % PM_BLK_SIZE;
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] >= 0)
                occupancy++;
            // if (testBit(&src.blk_list[blk_idx].myBitmap, bit_idx))
            //     occupancy++;
            right_index++;
            // bit_idx++;
            // if (bit_idx == PM_BLK_SIZE) {
            //     bit_idx = 0;
            //     blk_idx++;
            // }
        }
        density = (double)occupancy / (double)window_size;
        t_height = t_0 - (height * pma->delta_t);
        p_height = p_0 + (height * pma->delta_p);
        height++;
        
    } while ((density < p_height ||  
            density >= t_height) &&
           height < pma->tree_h);

    // cout << "height: " << height << " pma tree: " << pma->tree_h << endl;
    // print_pma(src);
    // print_vertex(src);
      /* Found a window within threshold. */
    if (height == 1)
        return;

    // cout << "test----" << endl;
    // int32_t tmpOccupacy = 0;
    // for (int32_t i = window_start; i < window_end; i++) {
    //     int32_t tmpBlkIdx = i / PM_BLK_SIZE;
    //     int32_t tmpBitIdx = i % PM_BLK_SIZE;
    //     if (src->blk_list[tmpBlkIdx].blk_ptr[tmpBitIdx] >= 0) {
    //         tmpOccupacy++;
    //     }
    // }
    // if (tmpOccupacy != occupancy)
    //     cout << "tmpOccupacy: " << tmpOccupacy << " occupancy: " << occupancy << endl;

    
    if (density >= p_height && density < t_height) {
        // pack (pma, window_start, window_end, occupancy);
        // cout << "test --- rebalance spread. occupancy: " << occupancy << endl;
        // print_pma(src);
        // print_vertex(src);
        spread (pma, src, window_start, window_end, occupancy, kind);
        // print_vertex(src);
    } else {
        // cout << "test --- rebalance resize. occupancy: " << occupancy << endl;
        resize_parallel (pma, src, thread_num, kind);
        // print_pma(src);
        // print_vertex(src);
  }
}

void resize_parallel(PMA* pma, Vertex *src, int thread_num, memkind_t kind)
{
    int32_t read_index = pma->arr_cap - 1;  
    int32_t readBitmapIdx = src->blk_list.size() - 1;
    // print_vertex(src);
    compute_capacity(pma);
    pma->tree_h = floor_lg (pma->seg_num) + 1;
    pma->delta_t = (t_0 - t_h) / pma->tree_h;
    pma->delta_p = (p_h - p_0) / pma->tree_h;


    int32_t blk_num = ceil((double)pma->arr_cap / PM_BLK_SIZE);  
    int32_t ex_blk_num = blk_num - src->blk_list.size();


    // cout << "ex_blk_num: " << ex_blk_num << endl;
    if (global_parallel_blk_vec[thread_num].size() < ex_blk_num) {
        generate_blk_batch_for_parallel(ex_blk_num * 2, thread_num, kind);
    }

    if (ex_blk_num > 0)
    {
        for (int32_t i=0; i < ex_blk_num; i++) {
            src->blk_list.push_back(global_parallel_blk_vec[thread_num].back());
            global_parallel_blk_vec[thread_num].pop_back();
        }
    } 

    // cout << "test -- 1" << endl;

    // print_vertex(src);
    // print_pma(src);

    int32_t capacity = pma->arr_cap;
    int32_t frequency = (capacity << 8) / pma->ele_num;
    int32_t write_index = (capacity << 8) - frequency;

    // cout << "read_idx: " << read_index << " write_index: " << (write_index >> 8) << endl;
    int32_t writeBitmapIdx = src->blk_list.size() - 1;
    
    while ((write_index >> 8) > read_index && read_index >= 0)
    {
        // cout << "read_idx: " << read_index << " write_index: " << (write_index >> 8) 
        //     << " readBitmapIdx: " << readBitmapIdx << " writeBitmapIdx: " << writeBitmapIdx << endl;
        int32_t write_off = (write_index >> 8) % PM_BLK_SIZE;
        int32_t read_off = read_index % PM_BLK_SIZE;
        
        // test
        // cout << "read_off: " << read_off << " write_off: " << write_off << endl;
        // src.blk_list[writeBitmapIdx].blk_ptr[write_off] = 1;
        // cout << "test ---  -1 " << endl;

        if (src->blk_list[readBitmapIdx].blk_ptr[read_off] >= 0) {
        // if (testBit(&src.blk_list[readBitmapIdx].myBitmap, read_off)) {
            // cout << "read_off: " << read_off << " write_off: " << write_off << endl;
            src->blk_list[writeBitmapIdx].blk_ptr[write_off] = src->blk_list[readBitmapIdx].blk_ptr[read_off];
            // cout << "test ---0 " << endl;
            // clearBit(&src.blk_list[readBitmapIdx].myBitmap, read_off);
            src->blk_list[readBitmapIdx].blk_ptr[read_off] = -1;
            // setBit(&src.blk_list[writeBitmapIdx].myBitmap, write_off);
            read_index--;
            write_index -= frequency;

            
            if ((read_index % PM_BLK_SIZE) > read_off) {
                readBitmapIdx --;
                // cout << "readBitmapIdx: " << readBitmapIdx << endl;
            }
            if (((write_index >> 8) % PM_BLK_SIZE) > write_off) {
                writeBitmapIdx--;
                // cout << "writeBitmapIdx: " << writeBitmapIdx << endl;
            }
        } else {
            // cout << "test ---1 " << endl;
            read_index--;
            if ((read_index % PM_BLK_SIZE) > read_off) {
                readBitmapIdx --;
                // cout << "readBitmapIdx: " << readBitmapIdx << endl;
            }
        }
    }
}

static void compute_capacity (PMA* pma) 
{
    /*，arr_cap*/
    pma->seg_len = ceil_lg (pma->ele_num);  /* Ideal segment size. */
    pma->seg_num = ceil_div (pma->ele_num, pma->seg_len);  /* Ideal number of segments. */
    /* The number of segments has to be a power of 2. */
    pma->seg_num = hyperceil (pma->seg_num);
    /* Update the segment size accordingly. */
    pma->seg_len = ceil_div (pma->ele_num, pma->seg_num);
    pma->arr_cap = pma->seg_len * pma->seg_num;
    /* Scale up as much as possible. */
    pma->arr_cap = max_sparseness * pma->arr_cap;
    pma->seg_len = max_sparseness * pma->seg_len;

    assert(pma->arr_cap > pma->ele_num);
}



// graph on NVM -- 20231116 wtf
void graph_init(int32_t max_vid, memkind_t kind)
{
	generate_blk_batch(max_vid * 3, kind);
	global_vectex_vec.resize(max_vid + 1);

    // cout << "test ---0" << endl;
	for (int32_t vidx = 0; vidx <= max_vid; vidx ++) {
        // cout << "vidx: " << vidx << " global_blk_vec size:" << global_blk_vec.size() << endl;
		Vertex* v = new Vertex(vidx);
		// v.id = vidx;
		// v.degree = 0;
		// v.blk_cnt = 1;
        // v.blk_list = (pm_blk *)malloc(sizeof(pm_blk) * v.blk_cnt);
		// v.blk_list[0] = global_blk_vec.back();
        v->blk_list.push_back(global_blk_vec.back());
        // cout << "initial test --- blk[0]: " << v.blk_list[0].blk_ptr[0] << endl;
		global_blk_vec.pop_back();
		v->pma = pma_creat(kind);
        v->edge_stamp_vec.resize(1);

        global_vectex_vec[vidx] = v;
	}
    cout << "graph initial complete!" << endl;
}

void graph_maintenance_nvm(memkind_t kind)
{
    cout << "Load graph on NVM baseline begin!" << endl;
	double start_t = get_current_time();
    // int64_t total_edge_num = global_edge_vec.size();
	for (int64_t edgeIdx = 0; edgeIdx < global_edge_num; edgeIdx++)
	{
        // if (edgeIdx == (total_edge_num / 2)) {
        //     start_t = get_current_time();
        // }

		// new_edge e = global_edge_vec[edgeIdx]; 
        new_edge e = global_edge_array[edgeIdx];  
		int32_t src = e.src;
        // int32_t src = 1;   
		int32_t des = e.des;
		Vertex* v = global_vectex_vec[src];
		int32_t index;
        // cout << "edgeIdx: " << edgeIdx << " src: " << src << " des: " << des << endl;
        
		// int32_t blk_idx = pma_find(global_vectex_vec[src].pma, global_vectex_vec[src], des, 0, (global_vectex_vec[src].blk_cnt-1), index);
		// index = blk_idx * PM_BLK_SIZE + index; 
        // cout << "blk_idx: " << blk_idx << " index: " << index << endl;
        index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1); // insert after index
        // cout << "test--" << endl;

        pma_insert(v->pma, v, des, index, kind);


        if ((edgeIdx && (edgeIdx % 10000000)) == 0) {
            // cout << "load " << (edgeIdx / 10000000) << " KW edges." << endl;
        }
    }
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    // print_graph();
    // print_pma(global_vectex_vec[1]);
	cout << "Load graph on NVM baseline time: " << end_t - start_t << endl;
    // print_graph();
}

void graph_maintenance_nvm_parallel_shuffle(memkind_t kind)
{
    cout << "Load graph on NVM baseline begin!" << endl;
	double start_t = get_current_time();
    // int64_t total_edge_num = global_edge_vec.size();
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for
	for (int64_t edgeIdx = 0; edgeIdx < global_edge_num; edgeIdx++)
	{
        // if (edgeIdx == (total_edge_num / 2)) {
        //     start_t = get_current_time();
        // }

		new_edge e = global_edge_vec[edgeIdx]; 
        // new_edge e = global_edge_array[edgeIdx];  
		int32_t src = e.src;
        // int32_t src = 1;  
		int32_t des = e.des;
		Vertex* v = global_vectex_vec[src];
        std::unique_lock<std::shared_timed_mutex> lock(v->mutex);
		int32_t index;
        // cout << "edgeIdx: " << edgeIdx << " src: " << src << " des: " << des << endl;

        index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1); // insert after index
        // cout << "test--" << endl;

        pma_insert(v->pma, v, des, index, kind);


        if ((edgeIdx && (edgeIdx % 10000000)) == 0) {
            // cout << "load " << (edgeIdx / 10000000) << " KW edges." << endl;
        }
    }
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    // print_graph();
    // print_pma(global_vectex_vec[1]);
	cout << "Load graph on NVM parallel shuffle time: " << end_t - start_t << endl;
    print_graph();
}

void graph_maintenance_nvm_parallel_shuffle_mix(memkind_t kind)
{
    cout << "Load graph on NVM parallel shuffle begin!" << endl;
	double start_t = get_current_time();
    // int64_t total_edge_num = global_edge_vec.size();
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for
	for (int64_t edgeIdx = 0; edgeIdx < global_edge_num; edgeIdx++)
	{
        // if (edgeIdx == (total_edge_num / 2)) {
        //     start_t = get_current_time();
        // }

		new_edge e = global_edge_vec[edgeIdx]; 
        // new_edge e = global_edge_array[edgeIdx];  
		int32_t src = e.src;
        // int32_t src = 1;  
		int32_t des = e.des;
		Vertex* v = global_vectex_vec[src];
        std::unique_lock<std::shared_timed_mutex> lock(v->mutex);
		int32_t index;
        // cout << "edgeIdx: " << edgeIdx << " src: " << src << " des: " << des << endl;

        index = binary_pma_find_mix(v->pma, v, des, 0, v->pma->arr_cap-1); // insert after index
        // cout << "test--" << endl;
        Vertex* v_src = global_vectex_vec[src];
        int status = e.statue;
        if (status < 0) { // delete edge
            continue;
            // int32_t blk_idx = index / PM_BLK_SIZE;  
            // int32_t bit_idx = index % PM_BLK_SIZE;
            // if (v_src->blk_list[blk_idx].blk_ptr[bit_idx] == des) {
            //     v_src->blk_list[blk_idx].blk_ptr[bit_idx] = -1;
            //     v_src->degree --;
            //     v_src->pma->ele_num--;
            // }
        } else {  // insert edge
            pma_insert(v_src->pma, v_src, des, index, kind);
        }
    }
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    // print_graph();
    // print_pma(global_vectex_vec[1]);
	cout << "Load graph on NVM parallel shuffle mix time: " << end_t - start_t << endl;
    print_graph();
}


void graph_maintenance_nvm_mv(memkind_t kind)
{
    // for edge insertion
    cout << "Load graph on NVM baseline begin!" << endl;
	double start_t = get_current_time();
    // int64_t total_edge_num = global_edge_vec.size();
	for (int64_t edgeIdx = 0; edgeIdx < global_edge_num; edgeIdx++)
	{
        // if (edgeIdx == (total_edge_num / 2)) {
        //     start_t = get_current_time();
        // }

		// new_edge e = global_edge_vec[edgeIdx]; 
        new_edge e = global_edge_array[edgeIdx];  
		int32_t src = e.src;
        // int32_t src = 1;  
		int32_t des = e.des;
		Vertex* v = global_vectex_vec[src];
		int32_t index;
        // cout << "edgeIdx: " << edgeIdx << " src: " << src << " des: " << des << endl;
        
		// int32_t blk_idx = pma_find(global_vectex_vec[src].pma, global_vectex_vec[src], des, 0, (global_vectex_vec[src].blk_cnt-1), index);
		// index = blk_idx * PM_BLK_SIZE + index; 
        // cout << "blk_idx: " << blk_idx << " index: " << index << endl;
        index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1); // insert after index
        // cout << "test--" << endl;

        pma_insert_mv(v->pma, v, des, index, kind);


        if ((edgeIdx && (edgeIdx % 10000000)) == 0) {
            // cout << "load " << (edgeIdx / 10000000) << " KW edges." << endl;
        }
    }
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);

	double end_t = get_current_time();
    // print_graph();
    // print_pma(global_vectex_vec[1]);
	cout << "Load graph on NVM baseline time: " << end_t - start_t << endl;
    print_graph();
}

void insert_edge_for_concurrent(new_edge e, memkind_t kind)
{

    writerCnt++;  

    int32_t src = e.src;
    // int32_t src = 1;  
    int32_t des = e.des;
    Vertex* v = global_vectex_vec[src];
    std::unique_lock<std::shared_timed_mutex> lock(v->mutex);
    int32_t index;

    index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1); // insert after index
    pma_insert_mv(v->pma, v, des, index, kind);
    

    // cout << "v1 next blk task id 1: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
    // cout << "v1 next blk task id 1.1: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
    // cout << "v1 next blk task id 1.2: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
    writerCnt--; 
    // cout << "insertion done!" << endl;
}

void graph_maintenance_nvm_parallel_baseline(Vertex *src, memkind_t kind)
{
    int32_t new_edges = src->edge_stamp_vec[0].size();
    for (int32_t edgeIdx = 0; edgeIdx < new_edges; edgeIdx++) {
        int32_t des = src->edge_stamp_vec[0][edgeIdx];
        int32_t index = binary_pma_find(src->pma, src, des, 0, src->pma->arr_cap-1);
        pma_insert(src->pma, src, des, index, kind);
    }
    src->edge_stamp_vec[0].clear();
    src->edge_stamp_vec[0].shrink_to_fit();
}

void graph_maintenance_nvm_insertion_parallel(Vertex *src, int64_t start_off, int64_t end_off, memkind_t kind)
{
    // int32_t new_edges = src->edge_stamp_vec[0].size();
    for (int64_t edgeIdx = start_off; edgeIdx < end_off; edgeIdx++) {
        int32_t des = global_edge_array[edgeIdx].des;
        int32_t index = binary_pma_find(src->pma, src, des, 0, src->pma->arr_cap-1);
        pma_insert(src->pma, src, des, index, kind);
    }
}

void graph_maintenance_nvm_mixed_parallel(Vertex *src, int64_t start_off, int64_t end_off, memkind_t kind)
{
    // int32_t new_edges = src->edge_stamp_vec[0].size();
    for (int64_t edgeIdx = start_off; edgeIdx < end_off; edgeIdx++) {
        int32_t des = global_edge_array[edgeIdx].des;
        int status = global_edge_array[edgeIdx].statue;
        int32_t index = binary_pma_find(src->pma, src, des, 0, src->pma->arr_cap-1);
        if (status < 0) { // delete edge
            int32_t blk_idx = index / PM_BLK_SIZE;  
            int32_t bit_idx = index % PM_BLK_SIZE;
            if (src->blk_list[blk_idx].blk_ptr[bit_idx] == des) {
                src->blk_list[blk_idx].blk_ptr[bit_idx] = -1;
                src->degree --;
                src->pma->ele_num--;
            } else {
                // cout << "error deletion: the edge didn't exist." << endl;
            }
        } else {  // insert edge
            pma_insert(src->pma, src, des, index, kind);
        }
    }
}

void graph_maintenance_nvm_parallel(Vertex *src, int thread_num, memkind_t kind)
{
    for (int32_t edgeIdx = 0; edgeIdx < src->edge_stamp_vec[0].size(); edgeIdx++) {
        int32_t des = src->edge_stamp_vec[0][edgeIdx];
        int32_t index = binary_pma_find(src->pma, src, des, 0, src->pma->arr_cap-1);
        pma_insert_parallel(src->pma, src, des, index, thread_num, kind);
    }
    src->edge_stamp_vec[0].clear();
}

void print_graph(void)
{
    double start_t = get_current_time();
    int64_t tottal_edges = 0;
    int64_t total_slots = 0;
    // cout << "v1 next blk task id 3: " << global_vectex_vec[1]->blk_list[0].next_blk->task_id << endl;
    for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        if (v->degree > 0) {
            // cout << vidx << ": ";
            for (int blk_idx = 0; blk_idx < v->blk_list.size(); blk_idx++) {
                // cout << "blk version: " << v->blk_list[blk_idx].task_id << endl;
                for (int ele_idx = 0; ele_idx < PM_BLK_SIZE; ele_idx++) {
                    if (v->blk_list[blk_idx].blk_ptr[ele_idx] >= 0) {
                        tottal_edges++;
                        // cout << v->blk_list[blk_idx].blk_ptr[ele_idx] << ", ";
                    } 
                    
                    total_slots++;
                }
            }
            // cout << endl;
        }
    }
    
    cout << "total number of edges: " << tottal_edges << " total slots: " << total_slots << endl;
    double end_t = get_current_time();
	cout << "Print graph time: " << end_t - start_t << endl;
}

void graph_maintenance_delete(memkind_t kind)
{
    cout << "Delete edges in batch!" << endl;
	double start_t = get_current_time();
    // int32_t total_edge_num = global_edge_vec.size();
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for
	for (int64_t edgeIdx = 0; edgeIdx < (global_edge_num/2); edgeIdx++)
	{
		// new_edge e = global_edge_vec[edgeIdx]; 
        new_edge e = global_edge_array[edgeIdx];  
		int32_t src = e.src;
        // int32_t src = 1;  
		int32_t des = e.des;
		Vertex* v = global_vectex_vec[src];
		int32_t index;
        index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
        int32_t blk_idx = index / PM_BLK_SIZE;  
        int32_t bit_idx = index % PM_BLK_SIZE;

        if (v->blk_list[blk_idx].blk_ptr[bit_idx] == des) {
            v->blk_list[blk_idx].blk_ptr[bit_idx] = -1;
            v->degree --;
            v->pma->ele_num--;
        } else {
            cout << "error deletion: the edge didn't exist." << endl;
        }

    }
    // cout << "delete over --" << endl;
    #pragma omp parallel for
    for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        if (v->degree == 0) { 
            #pragma omp critical
            {
                int32_t ori_blk_cnt = v->blk_list.size();
                for (int blk_idx = 1; blk_idx < ori_blk_cnt; blk_idx++) {
                    global_blk_vec.push_back(v->blk_list.back());
                    v->blk_list.pop_back();
                }
            }
            continue;
        }
        // cout << "test --0 degree: " << v->degree << " v->pma->arr_cap: " << v->pma->arr_cap  << endl;
        double density = (double)v->degree / (double)v->pma->arr_cap;
        if (v->blk_list.size() > 1 && (density < p_0)) { 
            PMA *pma = v->pma;
            int32_t read_index = pma->arr_cap - 1;  
            compute_capacity(pma);
            pma->tree_h = floor_lg (pma->seg_num) + 1;
            pma->delta_t = (t_0 - t_h) / pma->tree_h;
            pma->delta_p = (p_h - p_0) / pma->tree_h;
            
            // from:0 to:capacity
            int32_t capacity = pma->arr_cap;
            int32_t new_blk_num = capacity / PM_BLK_SIZE;  
            int32_t occupancy = v->degree;
            int32_t frequency = (capacity << 8) / occupancy; 
            int32_t write_index = (capacity << 8) - frequency;

            int32_t *tmp_win_arr = (int32_t*)malloc(capacity * sizeof(int32_t));

            while (read_index >= 0) { 
                int32_t blk_idx = read_index / PM_BLK_SIZE;
                int32_t bit_idx = read_index % PM_BLK_SIZE;

                if (v->blk_list[blk_idx].blk_ptr[bit_idx] > 0) {
                    tmp_win_arr[(write_index >> 8)] = v->blk_list[blk_idx].blk_ptr[bit_idx];
                    write_index -= frequency;
                    v->blk_list[blk_idx].blk_ptr[bit_idx] = -1;
                }
                read_index--;
            }
            int32_t ori_blk_cnt = v->blk_list.size();
            // cout << "test--1 capacity: " << capacity << " new_blk_num: " << new_blk_num << " ori_blk_cnt: " << ori_blk_cnt << endl;
            int32_t cpy_src = 0;
            int32_t cpy_des = 0;
            int32_t cpy_size = PM_BLK_SIZE;
            for (int32_t i = 0; i < new_blk_num; i++) {
                memcpy(&v->blk_list[i].blk_ptr[cpy_des], &tmp_win_arr[cpy_src], cpy_size * sizeof(int32_t));
                cpy_src += cpy_size;
                if (i == (new_blk_num - 2)) { 
                    cpy_size = (capacity-1) % PM_BLK_SIZE + 1; 
                }
            }
            // cout << "test--2" << endl;
            #pragma omp critical
            {
                for (int blk_idx = new_blk_num; blk_idx < ori_blk_cnt; blk_idx++) {
                    global_blk_vec.push_back(v->blk_list.back());
                    v->blk_list.pop_back();
                }
            }
            // cout << "test--3" << endl;
        }
    } 
    
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    // print_graph();
    // print_pma(global_vectex_vec[1]);
	cout << "Delete edges time: " << end_t - start_t << endl;
    // print_graph();
}

void delete_edge_for_concurrent(new_edge e, memkind_t kind)
{

    writerCnt++; 

    int32_t src = e.src;
    // int32_t src = 1;  
    int32_t des = e.des;
    Vertex* v = global_vectex_vec[src];
    std::unique_lock<std::shared_timed_mutex> lock(v->mutex); 

    int32_t index;

    index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1); // insert after index
    int32_t blk_idx = index / PM_BLK_SIZE;  
    int32_t bit_idx = index % PM_BLK_SIZE;


    if (v->blk_list[blk_idx].next_blk == NULL) {
        // cout << "blk task id: " << v->blk_list[blk_idx].task_id << " global_task_id: " << global_task_id << endl;
        if (v->blk_list[blk_idx].task_id < global_task_id) {
            // pm_blk blk = generate_blk(kind);  
            pm_blk* blk = new pm_blk(generate_blk(kind));
            blk->task_id = global_task_id;
            memcpy(blk->blk_ptr, v->blk_list[blk_idx].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
            v->blk_list[blk_idx].next_blk = blk;
        }
    } else {
        pm_blk* next_blk = v->blk_list[blk_idx].next_blk; 
        while (next_blk->next_blk != NULL) {
            next_blk = next_blk->next_blk;
        }
        if (next_blk->task_id < global_task_id) {  
            // pm_blk blk = generate_blk(kind);
            pm_blk* blk = new pm_blk(generate_blk(kind));
            blk->task_id = global_task_id;
            next_blk->next_blk = blk;
            memcpy(blk->blk_ptr, v->blk_list[blk_idx].blk_ptr, PM_BLK_SIZE * sizeof(int32_t));
        }
    }


    if (v->blk_list[blk_idx].blk_ptr[bit_idx] == des) {
        v->blk_list[blk_idx].blk_ptr[bit_idx] = -1;
        v->degree --;
        v->pma->ele_num--;
    } else {
        cout << "error deletion: the edge didn't exist." << endl;
    }

    writerCnt--; 
    // cout << "deletion done!" << endl;
}

void recycle_blk_for_del(vector<Vertex *> vset)
{

    #pragma omp parallel for
    for (int32_t vidx = 0; vidx < vset.size(); vidx++) {
        Vertex *v = vset[vidx];
        if (v->degree == 0) { 
            #pragma omp critical
            {
                int32_t ori_blk_cnt = v->blk_list.size();
                for (int blk_idx = 1; blk_idx < ori_blk_cnt; blk_idx++) {
                    global_blk_vec.push_back(v->blk_list.back());
                    v->blk_list.pop_back();
                }
            }
            continue;
        }

        double density = (double)v->degree / (double)v->pma->arr_cap;
        if (v->blk_list.size() > 1 && (density < p_0)) { 
            PMA *pma = v->pma;
            int32_t read_index = pma->arr_cap - 1;  
            compute_capacity(pma);
            pma->tree_h = floor_lg (pma->seg_num) + 1;
            pma->delta_t = (t_0 - t_h) / pma->tree_h;
            pma->delta_p = (p_h - p_0) / pma->tree_h;
            
            // from:0 to:capacity
            int32_t capacity = pma->arr_cap;
            int32_t new_blk_num = capacity / PM_BLK_SIZE;  
            int32_t occupancy = v->degree;
            int32_t frequency = (capacity << 8) / occupancy; 
            int32_t write_index = (capacity << 8) - frequency;

            int32_t *tmp_win_arr = (int32_t*)malloc(capacity * sizeof(int32_t));

            while (read_index >= 0) {
                int32_t blk_idx = read_index / PM_BLK_SIZE;
                int32_t bit_idx = read_index % PM_BLK_SIZE;

                if (v->blk_list[blk_idx].blk_ptr[bit_idx] > 0) {
                    tmp_win_arr[(write_index >> 8)] = v->blk_list[blk_idx].blk_ptr[bit_idx];
                    write_index -= frequency;
                    v->blk_list[blk_idx].blk_ptr[bit_idx] = -1;
                }
                read_index--;
            }

            
            int32_t cpy_src = 0;
            int32_t cpy_des = 0;
            int32_t cpy_size = PM_BLK_SIZE;
            for (int32_t i = 0; i < new_blk_num; i++) {
                memcpy(&v->blk_list[i].blk_ptr[cpy_des], &tmp_win_arr[cpy_src], cpy_size * sizeof(int32_t));
                cpy_src += cpy_size;
                if (i == (new_blk_num - 2)) {  
                    cpy_size = (capacity-1) % PM_BLK_SIZE + 1; 
                }
            }

            #pragma omp critical
            {
                int32_t ori_blk_cnt = v->blk_list.size();
                for (int blk_idx = new_blk_num; blk_idx < ori_blk_cnt; blk_idx++) {
                    global_blk_vec.push_back(v->blk_list.back());
                    v->blk_list.pop_back();
                }
            }
        }
    } 
}

void clean_snapshot(int32_t task_id, vector<int32_t> vset)
{

    #pragma omp parallel for 
    for (int32_t vidx = 0; vidx < vset.size(); vidx++) {  
        Vertex *v = global_vectex_vec[vset[vidx]];
        std::unique_lock<std::shared_timed_mutex> lock(v->mutex); 
        int32_t blk_cnt = v->blk_list.size();
        for (int32_t blk_idx = 0; blk_idx < blk_cnt; blk_idx++) {
            pm_blk *cur_blk = &(v->blk_list[blk_idx]);
            pm_blk *next_blk = cur_blk->next_blk;
            while (next_blk != NULL) {  
                if (next_blk->task_id <= task_id) {
                    cur_blk->next_blk = next_blk->next_blk;
                    // delete next_blk;
                    if (next_blk->task_id == task_id) {  
                        break;
                    }
                } else {  
                    break;
                }
                next_blk = next_blk->next_blk;
            }
        }
    }
}

void clean_snap_cur(void)
{
    while (1) {
        if (completeTaskVec[global_min_task] == 1) {
            vector<int32_t> vset;
            if (global_task_map.find(global_min_task) != global_task_map.end()) {
                for (const int32_t& element : global_task_map[global_min_task]) {
                    vset.push_back(element);
                }
                clean_snapshot(global_min_task, vset);
            }
            global_min_task++;
            if (global_min_task >= global_task_id) {  
                break;
            }
        }
    }
}

void print_vertex(Vertex *v)
{
    cout << "id: " << v->id << " { ";
    for (int32_t i=0; i < v->pma->arr_cap; i++) {
        int32_t blk_idx = i / PM_BLK_SIZE;
        int32_t bit_idx = i % PM_BLK_SIZE;
        // cout << v.blk_list[blk_idx].blk_ptr[bit_idx] << ", ";
        if (v->blk_list[blk_idx].blk_ptr[bit_idx] >= 0) {
            cout << v->blk_list[blk_idx].blk_ptr[bit_idx] << ", ";
        } else {
            // cout << -1 << ", ";
        }
    }
    cout << "}" << endl;
}

void print_vertex_mv(int32_t vid, Task task)
{
    Vertex *v = global_vectex_vec[vid];
    std::shared_lock<std::shared_timed_mutex> lock(v->mutex);  
    cout << "id: " << v->id << " { ";
      int32_t blk_cnt = v->blk_list.size();

        for (int32_t blk_idx = 0; blk_idx < blk_cnt; blk_idx++) {
            pm_blk *cur_blk = &(v->blk_list[blk_idx]);
            int32_t *blk_ptr = NULL;

            // blk_ptr = cur_blk->blk_ptr;
            // for (int32_t ele_idx = 0; ele_idx < PM_BLK_SIZE; ele_idx++) {
            //     // total_num_access++;
            //     if (blk_ptr[ele_idx] >= 0) {
            //         cout << blk_ptr[ele_idx] << ", ";
            //     }
            // }

            while (cur_blk->next_blk != NULL)
            {
                // cout << "test--- task id: " << cur_blk->task_id << " blk_cnt: " << blk_cnt
                //     << " cur_blk->task_id: " << cur_blk->task_id << endl;
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
                    cout << blk_ptr[ele_idx] << ", ";
                }
            }
        }
    cout << "}" << endl;
}

void print_pma(Vertex *v)
{
    cout << "seg_len: " << v->pma->seg_len << " ele_num: " << v->pma->ele_num << " tree_h: " << v->pma->tree_h <<
        " arr_cap: " << v->pma->arr_cap << " degree: " << v->degree << endl;
}

#endif
