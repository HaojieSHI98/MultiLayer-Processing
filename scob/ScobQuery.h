/*
 * HierKNN.h
 *
 *  Created on: 2016年9月8日
 *      Author: lostrong
 */

#ifndef SRC_HIERKNN_H_
#define SRC_HIERKNN_H_

#include<vector>
#include "../io/Input.h"
#include "../scob/ScobIndex.h"
#include "../RoadKNNUtilities.h"
using namespace std;


void prepare_environment(){
	srand(time(NULL));
	read_road_network();
	set_base_cell_num();
	read_scob_index();
}

void prepare_reverse_shortcut(int test_n, int tradeoff_level, int cut_level, int mode){
    for (int u = 0; u < test_n; u++)
        rev_shortcuts_arr[u].clear();

    for (int u = 0; u < test_n; u++) {

        vector<KNode> *nei;
        if (mode == 2) {
            if (node_order[u] >= cut_level) {
                nei = &(level_scs_arr[cut_level][u]);
                /*
                 * if node level >= cut_level, they are treated as with level $cut_level
                 */
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
                /*
                 * if at least tradeoff_level, then "level+up"
                 */
            else
                nei = &(shortcuts_inc_arr[u]);

        } else {
            if (node_order[u] >= cut_level) {
                continue;
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
            else
                nei = &(shortcuts_inc_arr[u]);

        }
        for (KNode knode: *nei) {
            KNode *node = new KNode(u, knode.dis);
            rev_shortcuts_arr[knode.id].push_back(*node);
        }

    }
}

void prepare_environment_with_scob_conf(int configurationId){
    prepare_environment();
    level_num = get_level_num_with_dataspace();
    cout<<"level_num: "<<level_num<<endl;
    vector<ScobConfiguration> configurations;
	if(network_name.compare("NY")==0)
		configurations = ScobConfiguration::generate_all_configurations(level_num);
	else
		configurations = ScobConfiguration::generate_all_update_favor_configurations(level_num);
    tradeoff_level = configurations[configurationId].tradeoff_level;
    cut_level = configurations[configurationId].cut_level;
    mode = configurations[configurationId].mode;
    prepare_reverse_shortcut(test_n, tradeoff_level, cut_level, mode);
}
vector<KNode> scob_KNN_query(int k, vector<KNode> *hier_local_knn, int query, long *dis, int *visited, DijkstraQueue *q) {

	int cnt_nodes = 0;
	int coor_x = cell_loc_x[query];
	int coor_y = cell_loc_y[query];

	vector<int> visitedVector;
	vector<KNode> result;
	dis[query] = 0;
	q->push(query, 0);
	visited[query] = 1;
	visitedVector.push_back(query);
	int naive_dijk = 0;
	if (cut_level == 0 && mode == 3) {
		naive_dijk = 1;
	}
	while (!q->empty()) {
		cnt_nodes++;
		int u = q->pop();
		if (naive_dijk) {
			//turns to naive Dijkstra
			if (result.size() >= k)
				break;
		} else if (result.size() >= k && dis[u] > result[k - 1].dis) {
			break;
		}
		int need_expand = 1;
		vector<KNode> &knodes = hier_local_knn[u];
		int size = knodes.size();
		if (naive_dijk) {
			KNode* node = new KNode(u, dis[u]);
			if (size)
				result.push_back(*node);
			if (result.size() >= k)
				break;
		} else {
			merge_k_list(result, knodes, k, dis[u]);
		}
		visited[u]=2;
		vector<KNode>* nei;

		if (mode == 2) {
			if (node_order[u] >= cut_level) {
				continue;
			}

			else if (node_order[u] >= tradeoff_level) {
				nei = &(shortcuts_inc_arr[u]);
			} else {
				nei = &(shortcuts_arr[u]);

			}


		} else {
			if (node_order[u] >= cut_level) {
				nei = &(level_scs_arr[cut_level][u]);
			}
			else if (node_order[u] >= tradeoff_level) {
				nei = &(shortcuts_inc_arr[u]);

			} else {
				nei = &(shortcuts_arr[u]);

			}


		}

		if (node_order[u]) {
			for (KNode knode : *nei) {
				int v = knode.id;
				if (visited[v] && dis[v] + knode.dis < dis[u]) {
					dis[u] = dis[v] + knode.dis;
					need_expand = 0;
					break;
				}
			}
		}

		if (need_expand) {
			for (KNode knode : *nei) {
				int v = knode.id;
				if (visited[v]) {

					if (knode.dis < dis[v] - dis[u]) {

						dis[v] = dis[u] + knode.dis;
						q->decrease_priority(v, dis[v]);
						visited[v] = 1;
					}
				} else {
					dis[v] = dis[u] + knode.dis;
					q->push(v, dis[v]);
					visited[v] = 1;
					visitedVector.push_back(v);
				}
			}

		}

	}

	for (int item : visitedVector)
		visited[item] = 0;
	q->clear();
	return result;
}

vector<KNode> scob_KNN_query_threshold(int k, vector<KNode> *hier_local_knn, int query, long *dis, int *visited, DijkstraQueue *q,
vector<int>& threshold, int query_id) {

    int cnt_nodes = 0;
    int coor_x = cell_loc_x[query];
    int coor_y = cell_loc_y[query];

    vector<int> visitedVector;
    vector<KNode> result;
    dis[query] = 0;
    q->push(query, 0);
    visited[query] = 1;
    visitedVector.push_back(query);
    int naive_dijk = 0;
    if (cut_level == 0 && mode == 3) {
        naive_dijk = 1;
    }
    while (!q->empty()) {
        cnt_nodes++;
        int u = q->pop();
        if (naive_dijk) {
            //turns to naive Dijkstra
            if (result.size() >= k)
                break;
        } else if (result.size() >= k && dis[u] > result[k - 1].dis) {
            break;
        }
        if(dis[u]>=threshold[query_id]) break;
        int need_expand = 1;
        vector<KNode> &knodes = hier_local_knn[u];
        int size = knodes.size();
        if (naive_dijk) {
            KNode* node = new KNode(u, dis[u]);
            if (size)
                result.push_back(*node);
            if (result.size() >= k)
                break;
        } else {
            merge_k_list(result, knodes, k, dis[u]);
        }
        visited[u]=2;
        vector<KNode>* nei;

        if (mode == 2) {
            if (node_order[u] >= cut_level) {
                continue;
            }

            else if (node_order[u] >= tradeoff_level) {
                nei = &(shortcuts_inc_arr[u]);
            } else {
                nei = &(shortcuts_arr[u]);

            }


        } else {
            if (node_order[u] >= cut_level) {
                nei = &(level_scs_arr[cut_level][u]);
            }
            else if (node_order[u] >= tradeoff_level) {
                nei = &(shortcuts_inc_arr[u]);

            } else {
                nei = &(shortcuts_arr[u]);

            }


        }

        if (node_order[u]) {
            for (KNode knode : *nei) {
                int v = knode.id;
                if (visited[v] && dis[v] + knode.dis < dis[u]) {
                    dis[u] = dis[v] + knode.dis;
                    need_expand = 0;
                    break;
                }
            }
        }

        if (need_expand) {
            for (KNode knode : *nei) {
                int v = knode.id;
                if (visited[v]) {

                    if (knode.dis < dis[v] - dis[u]) {

                        dis[v] = dis[u] + knode.dis;
                        q->decrease_priority(v, dis[v]);
                        visited[v] = 1;
                    }
                } else {
                    dis[v] = dis[u] + knode.dis;
                    q->push(v, dis[v]);
                    visited[v] = 1;
                    visitedVector.push_back(v);
                }
            }

        }

    }

    for (int item : visitedVector)
        visited[item] = 0;
    q->clear();
    return result;
}

vector<KNode> scob_KNN_query_edge(int k, vector<KNode> *hier_local_knn, int query, long *dis, int *visited,
								  DijkstraQueue *q) {
	int coor_x = cell_loc_x[query];
	int coor_y = cell_loc_y[query];
	vector<int> visitedVector;
	vector<KNode> result;
	dis[query] = 0;
	q->push(query, 0);
	visited[query] = 1;
	visitedVector.push_back(query);
	while (!q->empty()) {
		int u = q->pop();
		if (result.size() >= k && dis[u] > result[k - 1].dis) {
			break;
		}

		vector<KNode> knodes = hier_local_knn[u];


		int size = knodes.size();

		int t = 0;
		for (int i = 0; i < k && i < size; ++i) {
			KNode knode = knodes[i];
			int updated_dis = dis[u] + knode.dis;
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

			for (t2 = 0; t2 < result.size(); t2++) {
				if (result[t2].id == knode.id)
					break;
			}

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

		int need_expand = 1;
		vector<KNode> nei;

		if (mode == 1) {
			if (node_order[u] >= cut_level) {
				nei = level_scs_arr[cut_level][u];
			} else if (node_order[u] >= tradeoff_level)
				nei = shortcuts_arr[u];
			else
				nei = shortcuts_inc_arr[u];
		} else if (mode == 2) {
			if (node_order[u] >= cut_level) {
				continue;
			}

			else if (node_order[u] >= tradeoff_level)
				nei = shortcuts_inc_arr[u];
			else
				nei = shortcuts_arr[u];
		} else {
			if (node_order[u] >= cut_level) {
				nei = level_scs_arr[cut_level][u];
			} else if (node_order[u] >= tradeoff_level)
				nei = shortcuts_inc_arr[u];
			else
				nei = shortcuts_arr[u];
		}
		if (node_order[u] > 0) {
			for (KNode knode : nei) {
				int v = knode.id;
				if (visited[v] && dis[v] + knode.dis < dis[u]) {
					dis[u] = dis[v] + knode.dis;
					need_expand = 0;
					break;
				}
			}
		}
		if (need_expand) {
			for (KNode knode : nei) {
				int v = knode.id;
				int level_gap_for_v = level_gap[node_order[v]] * 2 + 1;
				int coor_v_x = cell_loc_x[v];
				int coor_v_y = cell_loc_y[v];
				int gap = max(abs(coor_x - coor_v_x), abs(coor_y - coor_v_y));
				if (gap > level_gap_for_v)
					continue;
				if (visited[v]) {
					if (knode.dis < dis[v] - dis[u]) {
						dis[v] = dis[u] + knode.dis;
						q->decrease_priority(v, dis[v]);
					}
				} else {
					dis[v] = dis[u] + knode.dis;
					q->push(v, dis[v]);
					visited[v] = 1;
					visitedVector.push_back(v);

				}

			}

		}

	}
	for (int i = 0; i < visitedVector.size(); i++) {
		visited[visitedVector[i]] = 0;
	}
	q->clear();
	return result;
}



#endif /* SRC_HIERKNN_H_ */
