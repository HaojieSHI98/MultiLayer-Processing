#ifndef UTILITY_CPP_
#define UTILITY_CPP_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
//#include <thread>
#include<unordered_map>
#include <limits.h>
//#include<omp.h>
#include <limits>
#include "dijk/DijkstraQueue.h"
#include "para/GlobalVariables.h"
#include "graph/Point.h"
#include "graph/Edge.h"
#include "graph/CellNode.h"
#include "io/Input.h"
#include "util/DataFunctions.h"


using namespace std;
extern int* x;
extern int* y;
extern vector<int>** cells;
extern vector<int>** circle_cells;
extern int* cell_loc;
extern int* circle_cell_loc[25];
extern int* regional_landmarks_all;
extern int eNext[2 * MAGIC_NUM + 2], eDis[2 * MAGIC_NUM + 1], eNode[2 * MAGIC_NUM + 1];
extern int *vStart;
extern int test_n;
extern int test_m;
extern int cell_num;
extern double largest_distance;
extern int x_min, y_min, x_max, y_max, x_len, y_len;
extern int* cell_loc_x;
extern int* cell_loc_y;

extern double EARTH_RADIUS;
extern double PI;






unsigned int my_rand() {
	return (rand() << 16) | (rand() & 0x0000ffff);
}



inline int test_in5(int i, int j, int u) {
	return abs(cell_loc[u] / cell_num - i) < 3 && abs(cell_loc[u] % cell_num - j) < 3;
}

inline int test_in3(int i, int j, int u) {
	return abs(cell_loc[u] / cell_num - i) < 2 && abs(cell_loc[u] % cell_num - j) < 2;
}

inline int test_in1(int i, int j, int u) {
	return abs(cell_loc[u] / cell_num - i) < 1 && abs(cell_loc[u] % cell_num - j) < 1;
}


void get_regional_landmarks(int h, int thread_num, int last_lev) {
	int *visited = new int[test_n + 1];
	memset(visited, 0, sizeof(int) * (test_n + 1));
	int *pre = new int[test_n + 1];
	memset(pre, 0, sizeof(int) * (test_n + 1));
	int *dis = new int[test_n + 1];
	for (int i = 1; i < test_n + 1; ++i)
		dis[i] = INT_MAX;
	int *all = new int[test_n + 1];
	int all_n = 0;
	int *it = new int[test_n + 1];
	int it_n = 0;
	size_t sum = 0;
	dis[0] = 0;
	DijkstraQueue* q = new DijkstraQueue(test_n);
	for (int i = h; i < cell_num; i += thread_num)
		for (int j = 0; j < cell_num; ++j) {
			vector<int> temp_landmarks1;
			vector<int> visited_list;
			int num = i * cell_num + j;
			for (int k = 0; k < cells[i][j].size(); ++k) {
				int s = cells[i][j][k];
				if (last_lev <= 0) {
					int b = 0;
					for (int x = vStart[s]; x != 0; x = eNext[x]) {
						if (cell_loc[eNode[x]] != num) {
							b = 1;
							break;
						}
					}
					if (!b)
						continue;
				} else {
					if (node_order[s] != last_lev)
						continue;
				}
				visited[s] = 1;
				dis[s] = 0;
				pre[s] = 0;
				q->clear();
				sum = circle_cells[i][j].size();
				q->push(s, 0);
				while (!q->empty() && sum != 0) {
					int u = q->pop();
//					inlev_order[u]++;
					all[all_n++] = u;
					if (test_in5(i, j, u)) {
						--sum;
					} else {
						it[it_n++] = u;
						continue;
					}
					for (int x = vStart[u]; x != 0; x = eNext[x]) {
						int v = eNode[x];
						if (!visited[v]) {
							pre[v] = u;
							visited[v] = 1;
							visited_list.push_back(v);
							dis[v] = dis[u] + eDis[x];
							q->push(v, dis[v]);
						} else {
							if (eDis[x] < dis[v] - dis[u]) {
								dis[v] = dis[u] + eDis[x];
								pre[v] = u;
								q->decrease_priority(v, dis[v]);
							}
						}
					}
				}
				while (!q->empty()) {
					int u = q->pop();
					all[all_n++] = u;
					it[it_n++] = u;
				}
				while (it_n > 0) {
					int u = it[--it_n];
					int flg = 1;
					while (pre[u] != 0) {
						if (!test_in3(i, j, u) && test_in3(i, j, pre[u])) {
							temp_landmarks1.push_back(pre[u]);
							if (test_in1(i, j, pre[u]) && !test_in3(i, j, u))
								temp_landmarks1.push_back(u);
							flg = 0;
						}

						int v = pre[u];
						pre[u] = 0;
						u = v;
					}
//					if (flg == 1)
//						cout << "error1" << endl;
//					if (!test_in1(i, j, u))
//						cout << "error2" << endl;

				}
				for (int ind : visited_list) {
					visited[ind] = 0;
					pre[ind] = 0;
					dis[ind] = INT_MAX;
				}
				all_n = 0;
				it_n = 0;

			}
			for (auto it = temp_landmarks1.cbegin(); it != temp_landmarks1.cend(); it++) {
				regional_landmarks_all[*it] = 1;
			}

		}

	delete q;
	delete[] it;
	delete[] all;
	delete[] dis;
	delete[] pre;
	delete[] visited;
}






void merge_k_list(vector<KNode>& result, vector<KNode>& mergelist, int k, int dis_u) {
	int t = 0;
	int size = mergelist.size();
	for (int i = 0; i < k && i < size; ++i) {
		KNode &knode = mergelist[i];
		int updated_dis = dis_u + knode.dis;
		if (result.size() >= k && updated_dis >= result[k - 1].dis) {
			break;
		}
		for (; t < result.size(); t++) {
			if (updated_dis < result[t].dis)
				break;
		}
		if (t == k) {
			break;
		}
		int t2 = 0;
		for (t2 = result.size() - 1; t2 >= 0; t2--) {
			if (result[t2].id == knode.id)
				break;
		}
		if (t2 == -1)
			t2 = result.size();
		if (t2 == result.size()) {
			if (result.size() < k) {
				KNode* c = new KNode(knode.id, updated_dis);
				result.push_back(*c);
			}
			for (int h = result.size() - 1; h > t; h--)
				result[h] = result[h - 1];
			if (t < k) {
				result[t].id = knode.id;
				result[t].dis = updated_dis;
			}

		} else {
			int t1 = max(t, t2);
			for (; t1 < result.size(); t1++) {
				if (result[t1].id == knode.id && updated_dis < result[t1].dis)
					break;
			}
			if (t1 < result.size()) {
				for (int h = t1; h > t; h--)
					result[h] = result[h - 1];
				if (t < k) {
					result[t].id = knode.id;
					result[t].dis = updated_dis;
				}
			}
		}

	}
}

void shuffle(vector<int> & values) {
    int val_sz = values.size();
    while (val_sz > 0) {
        int z1 = rand() % val_sz;
        int z2 = val_sz - 1;
        int tmp = values[z1];
        values[z1] = values[z2];
        values[z2] = tmp;
        val_sz--;
    }
}


void init_memory(){


	// init memory
	hier_local_knn_arr = new vector<KNode>[MAGIC_NUM];
	level_scs_arr = new vector<KNode> *[MAX_LEVELS];
	for (int i = 0; i < MAX_LEVELS; i++) {
		level_scs_arr[i] = new vector<KNode>[MAGIC_NUM];
		for (int j = 0; j < MAGIC_NUM; j++)
			level_scs_arr[i][j].clear();
	}


	shortcuts_inc_arr_local = new vector<KNode> *[MAX_LEVELS];
	for (int i = 0; i < MAX_LEVELS; i++) {
		shortcuts_inc_arr_local[i] = new vector<KNode>[MAGIC_NUM];
		for (int j = 0; j < MAGIC_NUM; j++) {
			shortcuts_inc_arr_local[i][j].clear();
		}
	}
	rev_shortcuts_arr = new vector<KNode>[MAGIC_NUM];
	rev_shortcuts_level_arr = new vector<KNode>[MAGIC_NUM];
}

void delete_memory(){
	delete[] hier_local_knn_arr;
	delete[] rev_shortcuts_arr;
	delete[] rev_shortcuts_level_arr;
//    cout<<"delete here1"<<endl;
	for (int i = 0; i < MAX_LEVELS; i++) {
		for (int j = 0; j < MAGIC_NUM; j++)
			level_scs_arr[i][j].clear();
		delete[] level_scs_arr[i];

	}
//    cout<<"delete here2"<<endl;
	delete[] level_scs_arr;

	for (int i = 0; i < MAX_LEVELS; i++) {
		for (int j = 0; j < MAGIC_NUM; j++)
			shortcuts_inc_arr_local[i][j].clear();
	}
	delete[] shortcuts_inc_arr_local;
//    cout<<"delete here3 "<<endl;
//    for (int i = 0; i < MAX_LEVELS; i++) {
//        for (int j = 0; j < MAGIC_NUM; j++)
//            shortcuts_inc_arr_local[i][j].clear();
//    }
//    delete[] shortcuts_inc_arr_local;
	cout<<"finish release memory..."<<endl;
	memory_deleted=1;

}

#endif /* UTILITY_CPP_ */
