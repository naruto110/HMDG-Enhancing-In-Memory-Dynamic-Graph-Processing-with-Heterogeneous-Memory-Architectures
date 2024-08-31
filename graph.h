#ifndef __GRAPH_H_
#define __GRAPH_H_

#include "global.h"

size_v findIdleBlk(vector<idle_blk*> &idle_blk_vec,size_v degree)
{
	size_v cap = degree + ceil((float)BUF_PRO * degree) + 3;  // +3 用来存储 id degree offset
	for (size_v blkIdx = 0; blkIdx < idle_blk_vec.size(); blkIdx++)
	{
		if (idle_blk_vec[blkIdx]->max_cap >= cap)
		{
			idle_blk_vec[blkIdx]->max_cap -= cap;  // 更新max_cap
			idle_blk_vec[blkIdx]->next_pos = idle_blk_vec[blkIdx]->next_pos - cap + 3;
			return idle_blk_vec[blkIdx]->blk_idx;
		}
	}
	return -1;
}

void toVertexV(const WorkerParams& params, char* line, PMEMoid root)
{
	// 数据格式：id |degree| neighbor-0 neighbor-1 ... 
	char* pch;
	pch = strtok(line, "\t");
	size_v id = atoi(pch);  // for larger id, use atol

	if (id > global_max_id)
	{
		global_max_id = id;
	}

	pch = strtok(NULL, " ");//skip neighbor_number
	size_v degree = atoi(pch);

	//cout << "ID: " << id << " degree: " << degree << endl;

	if (degree >= MAX_DEG_CAP) // 理论上一个blk装不下（包含buffer）
	{
		size_v numblk = ceil((float)degree / (MAX_LEN-3)); // 只考虑邻节点存储，不考虑buffer

		//cout << "numblk: " << numblk << endl;

		size_v nei_mod = degree % (MAX_LEN - 3);  // 最后一个blk存储的邻节点数量
		for (size_v i = 0; i < numblk; i++)
		{
			graph_blk* gblk = generate_new_blk(root);

			// 存储数据 前面blk存储，前面的blk不考虑buffer，最大存储 MAX_LEN - 3个neighbors
			gblk->vcnt = 1;
			gblk->buf[0] = id;
			
			gblk->buf[2] = MAX_LEN - 1;  // offset 从末位开始

			size_v neicnt = MAX_LEN - 3; // 对于单个顶点，单个blk最多能存邻节点数量
			size_v start_idx = MAX_LEN-1;  // 从后往前存储
			while ((neicnt--) && (pch = strtok(NULL, " "))) // 加载邻节点
			{
				size_v nb = atoi(pch);
				gblk->buf[start_idx] = nb;

				//cout << "index: " << start_idx << endl;

				start_idx--;
			}
			
			if (i < numblk - 1) // 除了最后一个blk，其他blk的next_blk保持-1
			{
				gblk->buf[1] = MAX_LEN - 3;  // 当前blk存储了MAX_LEN - 3个邻节点,即存满了
				gblk->next_blk = gblk->blk_id + 1;  // 确定知道还需要新建一个blk
			}
			else   // 最后一个blk
			{
				gblk->buf[1] = nei_mod;  //  最后一个blk未存满
			}
		}
	}
	else  // 一个blk能装下，优先从 idle_blk_vec 中选取blk进行存储
	{
		size_v blkIdx = findIdleBlk(idle_blk_vec, degree);
		if (blkIdx != -1)  // 找到合适的blk
		{
			graph_blk* gblk = graph_blk_vec[blkIdx];
			size_v start_idx = gblk->vcnt * 3;
			size_v preoff = gblk->buf[start_idx - 1];
			size_v predegree = gblk->buf[start_idx - 2];
			size_v prebuf = ceil((float)predegree * BUF_PRO);
			size_v offset = preoff - predegree - prebuf;  // 前一个顶点的起始offset - degree - buffer长度，即为当前顶点起始offset

			gblk->buf[start_idx] = id;
			gblk->buf[start_idx + 1] = degree;
			gblk->buf[start_idx + 2] = offset;

			while (pch = strtok(NULL, " ")) // 加载邻节点
			{
				size_v nb = atoi(pch);
				gblk->buf[offset] = nb;
				offset--;
			}
			gblk->vcnt++;

			size_v cur_cap = degree + ceil((float)degree * BUF_PRO);
			gblk->next_pos -= cur_cap;
		}
		else  // 新建一个新的blk
		{
			graph_blk* gblk = generate_new_blk(root);

			gblk->buf[0] = id;
			gblk->buf[1] = degree;
			size_v offset = MAX_LEN - 1;
			gblk->buf[2] = offset;  //  起始存储位置

			while (pch = strtok(NULL, " ")) // 加载邻节点
			{
				size_v nb = atoi(pch);
				gblk->buf[offset] = nb;
				offset--;
			}
			gblk->vcnt++;

			size_v cur_cap = degree + ceil((float)degree * BUF_PRO);
			gblk->next_pos -= cur_cap;

			size_v cap = ceil((float)degree * BUF_PRO) + degree + 3;
			size_v capmod = MAX_LEN - cap;

			idle_blk* iblk = (idle_blk*)malloc(sizeof(idle_blk));  // in DRAM

			iblk->blk_idx = gblk->blk_id;
			iblk->max_cap = capmod;
			idle_blk_vec.push_back(iblk);
		}
	}
}

void toEdge(const WorkerParams& params, char* line)
{
	// 数据格式：src_id des_id
	new_edge e;
	char* pch;
	pch = strtok(line, " ");
	size_v src = atoi(pch); 

	global_max_id = max(global_max_id, src);

	pch = strtok(NULL, " ");
	size_v des = atoi(pch);

	if (src == des)
		return;

	e.src = src; // 左侧连续，右侧离散
	e.des = des;
	// e.src = des;  // 右侧列连续，左侧离散
	// e.des = src;

	e.statue = 1; 
	global_edge_vec.push_back(e);
	global_edge_num++;

	e.src = des;
	e.des = src;
	global_edge_vec.push_back(e);
	global_edge_num++;
}

void toEdgeMix(const WorkerParams& params, char* line)
{
	// 数据格式：src_id des_id
	new_edge e;
	char* pch;
	pch = strtok(line, " ");
	size_v src = atoi(pch); 

	global_max_id = max(global_max_id, src);
	

	pch = strtok(NULL, " ");
	size_v des = atoi(pch);

	global_max_id = max(global_max_id, des);

	pch = strtok(NULL, " ");
	int status = atoi(pch);

	// int status = 0;
	// cout << status << " ";

	// cout << "status: " << status << " " << pch << " ";

	if (src == des)
		return;

	e.src = src; // 左侧连续，右侧离散
	e.des = des;
	e.statue = status;
	// e.src = des;  // 右侧列连续，左侧离散
	// e.des = src;

	// e.statue = 1; 
	global_edge_vec.push_back(e);
	global_edge_num++;

	// e.src = des;
	// e.des = src;
	// global_edge_vec.push_back(e);
	// global_edge_num++;
}

void load_graph(const WorkerParams& params, PMEMoid root)
{
	double start_t = get_current_time();
	string file_path;
	file_path = params.input_path;
	cout << " file: " << file_path << endl;

	ifstream ori_file(file_path.c_str());
	string line;
	int linecnt = 0;
	if (ori_file) // 有该文件  
	{
		while (getline(ori_file, line)) // line中不包括每行的换行符  
		{
			if (line.size() != 0)
			{
				linecnt++;
				//cout << "test: line---: " << linecnt << endl;
				toVertexV(params, const_cast<char*>(line.c_str()), root);
				// toEdge(params, const_cast<char*>(line.c_str()), root);
				// if (linecnt && linecnt % 10000000 == 0) {
				// 	cout << "load " << (linecnt / 10000000) << " edges" << endl;
				// }
			}
		}
		ori_file.close();
		//	cout << "worker: " << _my_rank << " load data over!!!" << endl;
	}
	else // 没有该文件  
	{
		cout << params.input_path.c_str() << endl;
		exit(-1);
	}

	double end_t = get_current_time();
	cout << "Load graph time: " << end_t - start_t << endl;
	// 构建顶点id与对应位置索引(很重要)
	//init(vertices);
}

void load_graph_dram(const WorkerParams& params)
{
	double start_t = get_current_time();
	string file_path;
	file_path = params.input_path;
	cout << " file: " << file_path << endl;

	ifstream ori_file(file_path.c_str());
	string line;
	int64_t linecnt = 0;
	if (ori_file) // 有该文件  
	{
		while (getline(ori_file, line)) // line中不包括每行的换行符  
		{
			if (line.size() != 0)
			{
				linecnt++;  // for test
				// if (linecnt == 26500)
				// {
				// 	return;
				// }
				// cout << "test: line---: " << linecnt << endl;
				// toVertexV(params, const_cast<char*>(line.c_str()), root);
				// toEdge(params, const_cast<char*>(line.c_str()));
				toEdgeMix(params, const_cast<char*>(line.c_str()));
			
			}
		}
		ori_file.close();
		cout << "total lines: " << linecnt << endl;
		//	cout << "worker: " << _my_rank << " load data over!!!" << endl;
	}
	else // 没有该文件  
	{
		cout << params.input_path.c_str() << endl;
		exit(-1);
	}

	double end_t = get_current_time();
	cout << "Load graph from disk time: " << end_t - start_t << endl;
	// 构建顶点id与对应位置索引(很重要)
	//init(vertices);
}

void print_csr(graph_blk* blk)
{
	cout << "graph print begin. Number of vertices: " << blk->vcnt << endl;
	int32_t startIdx = 0;
	int32_t offsetIdx = MAX_LEN;
	cout << "blk->buf: {";
	for (int32_t i=0; i<offsetIdx; i++)
	{
		cout << blk->buf[i] << ", ";
	}
	cout << "}" << endl;

	for (int32_t vidx=0; vidx < blk->vcnt; vidx++)
	{
		size_v vid = blk->buf[2*vidx];
		int32_t neiBeginIdx = blk->buf[2*vidx + 1];
		int32_t neiEndIdx = blk->buf[2*vidx + 3];

		// cout << "neiEndIdx: " << neiEndIdx-neiBeginIdx << endl;

		if (vidx == blk->vcnt-1)
		{
			neiEndIdx = blk->next_pos;
		}
		cout << "id: " << vid << " neiNum: " << neiBeginIdx - neiEndIdx << " neighbor list: {";
		for (int32_t neiIdx=neiBeginIdx; neiIdx > neiEndIdx; neiIdx--)
		{
			cout << blk->buf[neiIdx] << ", ";
		}
		cout << "}" << endl;
	}
	cout << "graph print complete." << endl;
}

void sorted_csr_pm(const WorkerParams& params, graph_blk* blk)  
{
	/*
		sorted定义：首先对比src大小，然后对比des大小
	*/
	double start_t = get_current_time();
	blk->blk_id = 0;
	blk->vcnt = 0;

	blk->next_pos = MAX_LEN - 1; // 从最后一位开始存储邻节点信息
	for (int edgeIdx = 0; edgeIdx < global_edge_vec.size(); edgeIdx++)
	{
		new_edge p = global_edge_vec[edgeIdx];
		size_v src = p.src;
		size_v des = p.des;

		// 在一个数组中实现CSR存储
		if (blk->vcnt == 0)  // 加入第一条边
		{
			int32_t startIdx = 2 * blk->vcnt;
			blk->buf[startIdx] = src;
			blk->buf[startIdx+1] = blk->next_pos;
			blk->buf[blk->next_pos] = des;

			startIdx += 2; // 每个顶点存储id和offset

			blk->vcnt += 1;
			blk->next_pos--; // 下次存储邻节点，--是因为邻节点存储从后往前
			continue;
		}

		// cout << "test---number of vertices: " << blk->vcnt << " id: " << blk->buf[0] << " " << blk->buf[MAX_LEN - 1] << endl;

		for (int32_t vidx = 0; vidx < blk->vcnt; vidx++)
		{
			/* 可替换为二分搜索src位置 */
			int32_t vid = blk->buf[2 * vidx];
			int32_t nextvid = blk->buf[2 * vidx + 2];
			int32_t nextoff = blk->buf[2 * vidx + 3];
			if (src == vid)  
			{
				/* 存在该src，只改变neighbor lish即可 */
				int32_t nextBeginIdx = blk->buf[2 * vidx + 3]; // 下一顶点的邻节点起始位置
				if (vidx < (blk->vcnt-1)) // 不是排在最后一个的src
				{
					int32_t tmpOff = blk->buf[2 * vidx + 1];  // 已存最后一个顶点的起始邻节点位置
					int32_t neiIdx = tmpOff;
					for (; neiIdx > nextBeginIdx; neiIdx--)
					{
						if (des < blk->buf[neiIdx])  // 从右往左，升序排列
						{
							//  所有后续元素后移一位
							for (int32_t j = blk->next_pos; j < neiIdx; j++)
							{
								blk->buf[j] = blk->buf[j+1];
							}
							blk->buf[neiIdx] = des;
							blk->next_pos--; // 下一个空位
							break;
						}
						// 需要插入当前顶点邻接表最后一个位置
						if ((neiIdx == nextBeginIdx+1) && (des > blk->buf[nextBeginIdx+1]))
						{
							for (int32_t j = blk->next_pos; j < nextBeginIdx; j++)
							{
								blk->buf[j] = blk->buf[j+1];
							}
							blk->buf[nextBeginIdx] = des;
							blk->next_pos--; // 下一个空位
						}
					}
					// 当前顶点之后所有顶点的Offset都需要-1，即左移一位
					for (int32_t remainvIdx = vidx+1; remainvIdx < blk->vcnt; remainvIdx++)
					{
						blk->buf[2*remainvIdx+1] -= 1;
					}
				}
				else // 当前src排在最后一个，直接在加入元素并排序即可. vidx = blk->vcnt-1
				{
					int32_t tmpOff = blk->buf[2 * vidx + 1];  // 已存最后一个顶点的起始邻节点位置
					int32_t neiIdx = tmpOff;
					for (; neiIdx > blk->next_pos; neiIdx--)
					{
						if (des < blk->buf[neiIdx])  // 从右往左，升序排列
						{
							//  所有后续元素后移一位
							for (int32_t j = blk->next_pos; j < neiIdx; j++)
							{
								blk->buf[j] = blk->buf[j+1];
							}
							blk->buf[neiIdx] = des;
							blk->next_pos--; // 下一个空位
							break;
						}

						if ((neiIdx == blk->next_pos+1) && (des > blk->buf[blk->next_pos+1]))
						{
							blk->buf[blk->next_pos] = des;
							blk->next_pos--; // 下一个空位
						}
					}
				}
				break;
			}
			else if (src>vid && vidx==(blk->vcnt-1))  // 无有效nextid
			{
				/* 当前不存在该src，新增src及其邻节点，只需注意src的排序，插在末尾 */
				blk->buf[2*blk->vcnt] = src;
				blk->buf[2*blk->vcnt+1] = blk->next_pos;

				blk->buf[blk->next_pos] = des;
				blk->next_pos--; // 下一个空位
				blk->vcnt += 1;
				break;
			}
			else if (src>vid && src<nextvid) // 存在有效nextid
			{
				/* 当前不存在该src，新增src及其邻节点，只需注意src的排序, 插在已有顶点中间 */
				// cout << "test--- vidx: " << vidx << " vcnt: " << blk->vcnt << endl;
				for (int32_t j = blk->vcnt; j>vidx; j--) // 移动顶点id和offset位置
				{
					blk->buf[2*j] = blk->buf[2*j-2];  // id
					blk->buf[2*j+1] = blk->buf[2*j-1]-1;  // offset，因为需要插入一个新值，需要再-1
				}

				blk->buf[2*vidx+2] = src;
				blk->buf[2*vidx+3] = nextoff;

				for (int32_t j = blk->next_pos; j<nextoff; j++)
				{
					blk->buf[j] = blk->buf[j+1];
				}
				blk->buf[nextoff] = des;
				blk->next_pos--; // 下一个空位
				blk->vcnt += 1;

				// print_csr(blk);
				break;
			}

		}
	}

	double end_t = get_current_time();
	// print graph
	// print_csr(blk);
	cout << "Load graph in CSR: " << end_t - start_t << endl;
}


void print_edges(void)
{
	cout << "{";
	for (int i = 0; i < global_edge_vec.size(); i++)
	{
		cout << "(" << global_edge_vec[i].src << ", " << global_edge_vec[i].des << ", " << global_edge_vec[i].statue << "), ";
	}
	cout << "}" << endl;
}

void generate_vertex_index_meta(vector<graph_blk*> graphBlkVec, vector<VIDX_PAIR> &vertexIdxMata)
{
	size_v total_blk = graphBlkVec.size();
	for (size_v blkIdx = 0; blkIdx < total_blk; blkIdx++)
	{
		size_v vertex_num = graphBlkVec[blkIdx]->vcnt;
		size_v bid = graphBlkVec[blkIdx]->blk_id;
		for (size_v vIdx = 0; vIdx < vertex_num; vIdx++)
		{
			size_v vid = graphBlkVec[blkIdx]->buf[vIdx * 3];
			VIDX_PAIR tmp_pair(blkIdx, vIdx);
			vertexIdxMata[vid] = tmp_pair;
		}
	}
}

size_v detect_buffer(graph_blk* gblk, size_v v_offset)
{
	/*
		验证一个顶点的buffer是否够的方法：
			1. 顶点offset - degree ？= next v offset
			2. 末端顶点需要特殊考虑
	*/
	size_v contributor_offset = -1;
	size_v hop = 1;
	size_v forward_flag = 1;
	size_v backward_flag = 1;
	size_v cur_offset = v_offset;

	while (true)
	{
		if (forward_flag == 1)
		{
			cur_offset = v_offset + hop;  // 拟从第cur_offset顶点处借buffer

			if (cur_offset == gblk->vcnt - 1)   // 最后一个顶点处是否有空余slot
			{
			//	if (gblk->buf[cur_offset * 3 + 2] - gblk->buf[cur_offset * 3 + 1] > cur_offset * 3 + 2)
				if (gblk->buf[cur_offset * 3 + 2] - gblk->buf[cur_offset * 3 + 1] > gblk->next_pos)
				{
					return cur_offset;
				}
			}
			else
			{
				if (cur_offset >= gblk->vcnt)  // 此处offset从0开始计数，到达vcnt时已经超出当前blk存储的总顶点数
				{
					forward_flag = -1;  // 不再往后找
				}
				else
				{
					if (gblk->buf[cur_offset * 3 + 2] - gblk->buf[cur_offset * 3 + 1] > gblk->buf[(cur_offset + 1) * 3 + 2])  // cur_offset有buffer
					{
						return cur_offset;
					}
				}
			}
		}

		if (backward_flag == 1)
		{
			cur_offset = v_offset - hop;
			if (cur_offset < 0)
			{
				backward_flag = -1;  // 不再往前找
			}
			else
			{
				if (gblk->buf[cur_offset * 3 + 2] - gblk->buf[cur_offset * 3 + 1] > gblk->buf[(cur_offset + 1) * 3 + 2])  // cur_offset有buffer
				{
					return cur_offset;
				}
			}
		}

		if (forward_flag == -1 && backward_flag == -1)
		{
			return -1;  // 没有空余buffer slot
		}
		hop++;
	}
}

void exchange_slot(PMEMobjpool* pop, graph_blk* gblk, size_v owner, size_v curv, size_v value)
{
	/*
		curv：当前需要插入邻节点的顶点（pos，gblk中的第几个顶点）
		owner：离curv最近的拥有空余buffer slot的顶点（pos）
		value：需要插入curv的邻节点id
	*/
	size_v stride; // slot左移还是右移
	size_v curoff = gblk->buf[curv * 3 + 2];
	size_v curdeg = gblk->buf[curv * 3 + 1];


	if (owner > curv) 
	{
		stride = 1;  // owner在当前顶点右侧，空余slot在当前顶点左侧
	}
	else
	{
		stride = -1; // owner在当前顶点左侧，空余slot在当前顶点右侧
	}

	if (stride == 1) // (邻节点列表在左侧)
	{
		while (owner != curv)  // move the slot
		{
			size_v owner_degree = gblk->buf[owner * 3 + 1];
			size_v owner_off = gblk->buf[owner * 3 + 2];
			gblk->buf[owner_off - owner_degree] = gblk->buf[owner_off];
			pmemobj_persist(pop, &(gblk->buf[owner_off - owner_degree]), sizeof(size_v));
			gblk->buf[owner * 3 + 2] = owner_off - 1;
			pmemobj_persist(pop, &(gblk->buf[owner * 3 + 2]), sizeof(size_v));

			owner -= stride;
		}
	}
	else
	{
		owner = owner + 1; // 左侧顶点（offset更大）从下一邻节点开始移动slot (邻节点列表在右侧)
		while (owner <= curv)
		{
			size_v owner_off = gblk->buf[owner * 3 + 2];
			size_v lastNeiPos = gblk->buf[owner * 3 + 2] - gblk->buf[owner * 3 + 1] + 1;  // owner最后一个邻节点存储的位置：offset - degree + 1

			gblk->buf[owner_off + 1] = gblk->buf[lastNeiPos];  // 最后一个邻节点复制到右侧一个空slot中
			pmemobj_persist(pop, &(gblk->buf[owner_off + 1]), sizeof(size_v));
			gblk->buf[owner * 3 + 2] = owner_off + 1;
			pmemobj_persist(pop, &(gblk->buf[owner * 3 + 2]), sizeof(size_v));

			owner -= stride;
		}
	}

	curoff = gblk->buf[curv * 3 + 2];
	gblk->buf[curoff - curdeg] = value;
	pmemobj_persist(pop, &(gblk->buf[curoff - curdeg]), sizeof(size_v));
	gblk->buf[curv * 3 + 1] = curdeg + 1;
	pmemobj_persist(pop, &(gblk->buf[curv * 3 + 1]), sizeof(size_v));
}

void update_vertex_index(vector<VIDX_PAIR>& vertexIdxMata, graph_blk* target_blk, size_v vid)
{
	size_v vertex_offset = target_blk->vcnt - 1;
	vertexIdxMata[vid].first = target_blk->blk_id;
	vertexIdxMata[vid].second = vertex_offset;
}

void remalloc_vertex(PMEMobjpool* pop, PMEMoid root, vector<VIDX_PAIR>& vertexIdxMata, vector<idle_blk*> &update_blk_vec, graph_blk* gblk, size_v vertex_index)
{
	/*
		将存储不下的顶点存储至其他blk中 (gblk中的第vertex_index个)
	*/
	size_v vid = gblk->buf[vertex_index * 3];
	size_v degree = gblk->buf[vertex_index * 3 + 1];
	size_v ori_offset = gblk->buf[vertex_index * 3 + 2];

	size_v blkIdx = findIdleBlk(update_blk_vec, degree);  // 此处找buffer需要重新更新最后一个顶点满buffer状态
//	size_v blkIdx = findIdleBlkInsert(idle_blk_map, degree);
	if (blkIdx != -1)  // 找到合适的blk
	{
		graph_blk* target_blk = graph_blk_vec[blkIdx];
		size_v start_idx = target_blk->vcnt * 3;  // 顶点数作为索引位置
		size_v preoff = target_blk->buf[start_idx - 1];
		size_v predegree = target_blk->buf[start_idx - 2];

	//	size_v prebuf = ceil((float)predegree * BUF_PRO);  // 此处buffer 不好确定
	//	size_v offset = preoff - predegree - prebuf;  // 前一个顶点(最后一个顶点)的起始preoff - degree - buffer长度，即为当前顶点起始offset
		size_v offset = target_blk->next_pos;  // 其实填充邻节点的位置

		target_blk->buf[start_idx] = vid;
		pmemobj_persist(pop, &(target_blk->buf[start_idx]), sizeof(size_v));
		target_blk->buf[start_idx + 1] = degree;
		pmemobj_persist(pop, &(target_blk->buf[start_idx + 1]), sizeof(size_v));
		target_blk->buf[start_idx + 2] = offset;  // 在该blk的offset处开始存储邻节点
		pmemobj_persist(pop, &(target_blk->buf[start_idx + 2]), sizeof(size_v));

		for (size_v neiIdx = 0; neiIdx < degree; neiIdx++)
		{
			target_blk->buf[offset - neiIdx] = gblk->buf[ori_offset - neiIdx];
			pmemobj_persist(pop, &(target_blk->buf[offset - neiIdx]), sizeof(size_v));
		}

		target_blk->vcnt++;
		pmemobj_persist(pop, &(target_blk->vcnt), sizeof(size_v));

		size_v tmp_cap = degree + ceil((float)BUF_PRO * degree);
		target_blk->next_pos -= tmp_cap;  // 这个不需要persistent

		update_vertex_index(vertexIdxMata, target_blk, vid);  // 更新顶点索引

		
		//print_array(target_blk->buf, MAX_LEN);
	}
	else  // 新建一个新的blk
	{
		graph_blk* target_blk = generate_new_blk(root);

		target_blk->buf[0] = vid;
		pmemobj_persist(pop, &(target_blk->buf[0]), sizeof(size_v));
		target_blk->buf[1] = degree;
		pmemobj_persist(pop, &(target_blk->buf[1]), sizeof(size_v));

		size_v offset = MAX_LEN - 1;
		target_blk->buf[2] = offset;  //  起始存储位置
		pmemobj_persist(pop, &(target_blk->buf[2]), sizeof(size_v));

		for (size_v neiIdx = 0; neiIdx < degree; neiIdx++)
		{
			target_blk->buf[offset - neiIdx] = gblk->buf[ori_offset - neiIdx];
			pmemobj_persist(pop, &(target_blk->buf[offset - neiIdx]), sizeof(size_v));
		}

		target_blk->vcnt++;
		pmemobj_persist(pop, &(target_blk->vcnt), sizeof(size_v));
		

		size_v cap = ceil((float)degree * BUF_PRO) + degree + 3;
		size_v capmod = MAX_LEN - cap;

		target_blk->next_pos = MAX_LEN - 1 - (cap - 3);

		//idle_blk_map[target_blk->blk_id] = capmod;

		idle_blk* iblk = (idle_blk*)malloc(sizeof(idle_blk));  // in DRAM

		iblk->blk_idx = target_blk->blk_id;
		iblk->max_cap = capmod; // 可能为负数
		iblk->next_pos = MAX_LEN - 1 - (cap - 3);  //  下一个可存储邻节点的位置
		update_blk_vec.push_back(iblk);

		//print_array(target_blk->buf, MAX_LEN);
		//cout << target_blk->next_pos << endl;

		update_vertex_index(vertexIdxMata, target_blk, vid);  // 更新顶点索引
	}
}

void update_blk(PMEMobjpool* pop, graph_blk* gblk, size_v max_vidx, size_v buffer_cnt)
{
	/*
		对于前max_vidx个顶点进行buffer重分配
		从最后一个顶点开始移动
		移动完成后更新当前blk存储的顶点数
	*/
	gblk->vcnt = max_vidx + 1;
	pmemobj_persist(pop, &(gblk->vcnt), sizeof(size_v));

	for (size_v vidx = max_vidx; vidx >= 0; vidx--)
	{
		size_v degree = gblk->buf[vidx * 3 + 1];
		size_v offset = gblk->buf[vidx * 3 + 2];
		size_v lastNeiPos = offset - degree + 1;  // 最后一个邻节点存储的位置
		size_v tmpBufSize = ceil((float)degree * BUF_PRO);  // 当前顶点所需要的buffer size
		buffer_cnt -= tmpBufSize;

		size_v new_offset = offset - buffer_cnt;  // 空出buffer给剩余顶点 （最后一个顶点的邻节点整体左移buffer_cnt个slot）
		// 不一定需要全部移动，因为邻节点数量可能大于剩余buffer长度，此情况只需移动部分邻节点
		if (new_offset > lastNeiPos)
		{
			for (size_v neiIdx = 0; neiIdx < buffer_cnt; neiIdx++)
			{
				gblk->buf[lastNeiPos - 1 - neiIdx] = gblk->buf[offset - neiIdx];  // 只需要把前重合部分移动只buffer区即可，避免数据覆盖
				pmemobj_persist(pop, &(gblk->buf[lastNeiPos - 1 - neiIdx]), sizeof(size_v));
			}
		}
		else   // 所有邻节点都需要移动
		{
			for (size_v neiIdx = 0; neiIdx < degree; neiIdx++)
			{
				gblk->buf[new_offset - neiIdx] = gblk->buf[offset - neiIdx];
				pmemobj_persist(pop, &(gblk->buf[new_offset - neiIdx]), sizeof(size_v));
			}
		}
		gblk->buf[vidx * 3 + 2] = new_offset;  // 所有邻节点移动完毕 更新offset
		pmemobj_persist(pop, &(gblk->buf[vidx * 3 + 2]), sizeof(size_v));

		if (vidx == max_vidx)
		{
			gblk->next_pos = new_offset - degree - tmpBufSize;  // 最后一个顶点确定位置后即可确定当前blk的next_pos
		}
	}

	//print_array(gblk->buf, MAX_LEN);
/*	gblk->vcnt = max_vidx + 1;
	pmemobj_persist(pop, &(gblk->vcnt), sizeof(size_v));*/   // 移动完成后更新当前blk存储的顶点数
}

void generate_and_split(PMEMobjpool* pop, PMEMoid root, vector<VIDX_PAIR>& vertexIdxMata, vector<idle_blk*>& update_blk_vec, graph_blk* gblk, size_v src, size_v des)
{
	/*
		step 0: 确定当前gblk能存储的最多顶点数
		step 1: 将其他顶点分配到其他有空闲的blk中
		step 2: 无空余blk则申请新的blk，并将顶点分配至新blk中
	*/
	size_v total_occupation = 0;
	size_v buffer_cnt = 0;
	size_v max_vidx = 0;
	for (size_v vidx = 0; vidx < gblk->vcnt; vidx++)
	{
		size_v degree = gblk->buf[vidx * 3 + 1];
		size_v occupation = 3 + degree + ceil((float)degree * BUF_PRO);
		total_occupation += occupation;
		if (total_occupation > MAX_LEN)
		{
			max_vidx = vidx - 1;  // 只能满足到前一个顶点的存储需求
			break;
		}
		buffer_cnt += ceil((float)degree * BUF_PRO);
	}

	//cout << "maximum occupation: " << max_vidx << endl;

	for (size_v vidx = max_vidx + 1; vidx < gblk->vcnt; vidx++)
	{
		remalloc_vertex(pop, root, vertexIdxMata, update_blk_vec, gblk, vidx);
	}

	//cout << "remalloc_vertex complete !!!" << endl;

	// 当所有顶点重新分配完成后，更新原gblk中的顶点数，并重新分配空余slot(按需分配),即第0-max_vidx的顶点进行buffer重分配（根据BUF_PRO）
	update_blk(pop, gblk, max_vidx, buffer_cnt);

	//print_array(gblk->buf, MAX_LEN);
	//cout << "update_blk complete !!!" << endl;

	// 直接插入新邻节点，此时自己buffer充足，直接插入即可
	VIDX_PAIR srcp = vertexIdxMata[src];
	size_v blk_id = srcp.first;
	size_v v_offset = srcp.second;  // 存储的第几个顶点
	graph_blk* gblk_new = graph_blk_vec[blk_id];
	size_v vidx = v_offset * 3;
	size_v degree = gblk_new->buf[vidx + 1];
	size_v nei_offset = gblk_new->buf[vidx + 2];

	gblk_new->buf[nei_offset - degree] = des;
	pmemobj_persist(pop, &(gblk_new->buf[nei_offset - degree]), sizeof(size_v));
	gblk->buf[vidx + 1] += 1;  // degree加1
	pmemobj_persist(pop, &(gblk->buf[vidx + 1]), sizeof(size_v));

}

void super_vertex_insertion(PMEMobjpool* pop, PMEMoid root, graph_blk* gblk, size_v value)
{
	/*
	* gblk：存储需要插入新邻节点的顶点所在的blk
	* value：待插入的值
		1. 根据当前blk所存储的邻节点数量判定是否从当前blk插入新邻节点
		2. 当前blk存满了，根据next_blk从下一blk中寻找空余slot
		3. 查找到最后一个blk，即next_blk为-1，弱也存满了则新建一个新blk，并讲当前blk的next_blk修改为对应新blk号
	*/
	size_v next_blk;
	size_v degree;
	graph_blk* cur_blk = gblk;

	while (cur_blk->next_blk)
	{
		degree = cur_blk->buf[1];
		if (degree == MAX_LEN - 3)  // 存满了
		{
			next_blk = cur_blk->next_blk;
			cur_blk = graph_blk_vec[next_blk];
		}
		else
		{
			break;  // 当前cur_blk存在空余buffer
		}
	}

	if (cur_blk->next_blk == -1)  // 最后一个blk，进一步判断是否有空余slot
	{
		degree = cur_blk->buf[1];
		if (degree == MAX_LEN - 3)   // 存满了，需要新建新的blk
		{
			graph_blk* new_blk = generate_new_blk(root);
			new_blk->buf[0] = cur_blk->buf[0];
			new_blk->buf[1] = 1;
			new_blk->buf[2] = MAX_LEN - 1;  // 末位
			new_blk->buf[MAX_LEN - 1] = value;
			pmemobj_persist(pop, new_blk, sizeof(graph_blk));
			
			cur_blk->next_blk = new_blk->blk_id;
			pmemobj_persist(pop, &(cur_blk->next_blk), sizeof(size_v));
		}
		else
		{
			size_v offset = cur_blk->buf[2];
			cur_blk->buf[offset - degree] = value;
			pmemobj_persist(pop, &(cur_blk->buf[offset - degree]), sizeof(size_v));
			cur_blk->buf[1] += 1;
			pmemobj_persist(pop, &(cur_blk->buf[1]), sizeof(size_v));
		}
	}
	else   // 必然有空余slot
	{
		size_v offset = cur_blk->buf[2];
		cur_blk->buf[offset - degree] = value;
		pmemobj_persist(pop, &(cur_blk->buf[offset - degree]), sizeof(size_v));
		cur_blk->buf[1] += 1;
		pmemobj_persist(pop, &(cur_blk->buf[1]), sizeof(size_v));
	}
}

void insert_edge(PMEMobjpool* pop, PMEMoid root, vector<VIDX_PAIR>& vertexIdxMata, vector<idle_blk*>& update_blk_vec, size_v src, size_v des)
{
	/*
		* pop：内存池，用于持久存储数据
		* root：用于新建blk
		* vertexIdxMata：根据顶点id查询顶点所在blk及偏移量v_offset
		* update_blk_vec: 用于增边扩充的新blk集合
		* src/des: 将des插入src的邻节点列表中
	
		* insert des into src
			case 1: enough buffer
			case 2: borrow from neighbors (maybe n-hop neighbors)
			case 3: split block (re-malloc buffer according to BUF_PRO)
			case 4: 针对单顶点无法存储在单个blk的情况 （待完成）
	*/
	VIDX_PAIR srcp = vertexIdxMata[src];
	size_v blk_id = srcp.first;
	size_v v_offset = srcp.second;  // 存储的第几个顶点
	graph_blk* gblk = graph_blk_vec[blk_id];
	size_v vidx = v_offset * 3;
	size_v degree = gblk->buf[vidx + 1];
	size_v nei_offset = gblk->buf[vidx + 2];
	size_v buffer_size = 0;
	size_v next_offset = gblk->buf[vidx + 5];

	// 判断是否是超级顶点（多个blk存储）
	if (gblk->next_blk != -1)
	{
		cout << "case 0: super vertex!" << endl;

		super_vertex_insertion(pop, root, gblk, des);
		return;
	}

	if (v_offset + 1 < gblk->vcnt) // 不是最后一个顶点
	{
		buffer_size = nei_offset - next_offset - degree;
	}
	else
	{
	//	buffer_size = nei_offset - degree - (v_offset + 2);  // 最后一个顶点插入了新邻节点
		buffer_size = nei_offset - degree - gblk->next_pos;  //  最后一个顶点的buffer
	}

	if (buffer_size > 0)  // case 1
	{
		cout << "case 1: have enough buffer slots!" << endl;

		size_v insert_pos = nei_offset - degree;  

		gblk->buf[insert_pos] = des;
		pmemobj_persist(pop, &(gblk->buf[insert_pos]), sizeof(size_v));

		gblk->buf[vidx + 1] += 1;  // degree加1
		pmemobj_persist(pop, &(gblk->buf[vidx + 1]), sizeof(size_v));
	}
	else  // buffer不够了
	{
		// case 2: borrow from neighbors (maybe n-hop neighbors)

		size_v owner_pos = detect_buffer(gblk, v_offset);
		cout << "owner_pos: " << owner_pos << endl;

		if (owner_pos != -1)  // 找到空余slot
		{
			cout << "case 2: borrow slots from neighbors!" << endl;

			exchange_slot(pop, gblk, owner_pos, v_offset, des);
		}
		else  // case 3: require new blk -> split the original blk
		{
			cout << "case 3: generate new block!" << endl;

			generate_and_split(pop, root, vertexIdxMata, update_blk_vec, gblk, src, des);
		}
	}

}

void delete_edge(PMEMobjpool* pop, size_v src, size_v des, vector<VIDX_PAIR> vertexIdxMata)
{
	/*
		delete des from src
		将要删除的邻节点与最后一个邻节点更换位置，并修改offset
		（注意区别对待超级顶点）
	*/
	VIDX_PAIR srcp = vertexIdxMata[src];
	size_v blk_id = srcp.first;
	size_v v_offset = srcp.second;  // 存储的第几个顶点
	graph_blk* gblk = graph_blk_vec[blk_id];
	size_v vidx = v_offset * 3;
	size_v degree = gblk->buf[vidx + 1];
	size_v nei_offset = gblk->buf[vidx + 2];

	// 处理超级顶点：在被删除邻节点所在blk进行删除即可，需更改当前blk的degree值
	if (gblk->next_blk != -1)
	{
		graph_blk* cur_blk = gblk;
		degree = cur_blk->buf[1];
		nei_offset = cur_blk->buf[2];
		size_v break_flag = 0;
		size_v next_blk = cur_blk->next_blk;
		while (break_flag == 0)
		{
			for (size_v neiIdx = 0; neiIdx < degree; neiIdx++)
			{
				if (cur_blk->buf[nei_offset - neiIdx] == des)
				{
					cur_blk->buf[nei_offset - neiIdx] = cur_blk->buf[nei_offset - degree + 1]; // 将末尾邻节点复制至需要删除的邻节点处
					pmemobj_persist(pop, &(cur_blk->buf[nei_offset - neiIdx]), sizeof(size_v));
					cur_blk->buf[vidx + 1] -= 1;  // 当前blk degree减1
					pmemobj_persist(pop, &(cur_blk->buf[vidx + 1]), sizeof(size_v));  // recovery是时候需要额外判定是否存在与最后一个邻节点id相同的邻节点，如果有则说明此处flush时crash了
					
					break_flag = 1;
					break;
				}
			}

			if (break_flag == 0)
			{
				cur_blk = graph_blk_vec[next_blk];
				degree = cur_blk->buf[1];
				nei_offset = cur_blk->buf[2];
				next_blk = cur_blk->next_blk;
			}
		}

		return;
	}


	// 处理普通顶点，即单个blk存储了多个顶点
	for (size_v neiIdx = 0; neiIdx < degree; neiIdx++)
	{
		if (gblk->buf[nei_offset - neiIdx] == des)
		{
			gblk->buf[nei_offset - neiIdx] = gblk->buf[nei_offset - degree + 1]; // 将末尾邻节点复制至需要删除的邻节点处
			pmemobj_persist(pop, &(gblk->buf[nei_offset - neiIdx]), sizeof(size_v));
			gblk->buf[vidx + 1] -= 1;
			pmemobj_persist(pop, &(gblk->buf[vidx + 1]), sizeof(size_v));  // recovery是时候需要额外判定是否存在与最后一个邻节点id相同的邻节点，如果有则说明此处flush时crash了
			return;
		}
	}
}

void insert_single_edge(PMEMobjpool* pop, PMEMoid root, vector<VIDX_PAIR>& vertexIdxMata, vector<idle_blk*>& update_blk_vec, newEdge edge)
{
	/*
		case 1: enough buffer
		case 2: borrow from neighbors
		case 3: split block (re-malloc buffer according to BUF_PRO)
	*/
	size_v srcv = edge.src;
	size_v desv = edge.des;

	insert_edge(pop, root, vertexIdxMata, update_blk_vec, srcv, desv);
}

void insert_batch_edges(void)
{

}

void delete_single_edge(PMEMobjpool* pop, vector<VIDX_PAIR> vertexIdxMata, newEdge edge)
{
	size_v srcv = edge.src;
	size_v desv = edge.des;

	delete_edge(pop, srcv, desv, vertexIdxMata);
}

void graph_recovery(void)
{
	
}

void kcore_recovery(void)
{

}




#endif
