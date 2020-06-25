//
// Created by lostrong on 17/7/10.
//

#ifndef SOB_TESTSHORTESTPATH_H
#define SOB_TESTSHORTESTPATH_H
#include "../dijk/DijkstraQueue.h"
#include "../para/GlobalVariables.h"
#include "../graph/CellNode.h"
#include "../util/MemoryAllocation.h"
int TestShortestPath(int s, int t, int *node_order, unordered_map<int, vector<KNode> > &shortcuts, mem_struct &mems,
                     int thresh_dist) {
    int min_length = INT_MAX;
    DijkstraQueue *q1 = mems.q;
    DijkstraQueue *q2 = mems.q1;
    long *dis1 = mems.dist;
    long *dis2 = mems.dist1;
    int *visited1 = mems.visited;
    int *visited2 = mems.visited1;
    q1->clear();
    q2->clear();
    vector<int> visitedVector1;
    vector<int> visitedVector2;
    int valid1 = 1;
    int valid2 = 1;
    dis1[s] = 0;
    dis2[t] = 0;
    q1->push(s, 0);
    q2->push(t, 0);
    visited1[s] = 1;
    visited2[t] = 1;
    visitedVector1.push_back(s);
    visitedVector2.push_back(t);
//	cout<<"there"<<endl;
    while ((!q1->empty() && valid1) || (!q2->empty() && valid2)) {
        if (!q1->empty() && valid1) {
            int u1 = q1->pop();
            if (dis1[u1] > min_length || dis1[u1] > thresh_dist) {
                valid1 = 0;
            }
            visited1[u1] = 2;

            int need_expand = 1;
            vector<KNode> nei = shortcuts_arr[u1];
            for (KNode knode : nei) {
                int v = knode.id;
                if (visited1[v] && dis1[v] + knode.dis < dis1[u1]) {
                    dis1[u1] = dis1[v] + knode.dis;
                    need_expand = 0;
                    break;
                }
            }
            if (visited2[u1]) {
                if (min_length > dis1[u1] + dis2[u1])
                    min_length = dis1[u1] + dis2[u1];

            }

            if (need_expand) {
                for (KNode knode : nei) {
                    int v = knode.id;
//					if (node_order[v] < node_order[u1]
//							|| (node_order[v] == node_order[u1]
//									&& inlev_order[v] < inlev_order[u1])
//							|| (node_order[v] == node_order[u1]
//									&& inlev_order[v] == inlev_order[u1]
//									&& v < u1)) {
                    if (node_order[v] < node_order[u1]) {
                        continue;
                    }
                    if (visited1[v]) {
                        if (knode.dis < dis1[v] - dis1[u1]) {
                            dis1[v] = dis1[u1] + knode.dis;
                            q1->decrease_priority(v, dis1[v]);
                        }
                    } else {
                        dis1[v] = dis1[u1] + knode.dis;
                        q1->push(v, dis1[v]);
                        visited1[v] = 1;
                        visitedVector1.push_back(v);

                    }

                }
            }
        }

        if (!q2->empty() && valid2) {
            int u2 = q2->pop();
            if (dis2[u2] > min_length || dis2[u2] > thresh_dist) {
                valid2 = 0;
            }
            visited2[u2] = 2;

            int need_expand = 1;
            vector<KNode> nei = shortcuts_inc_arr[u2];
            for (KNode knode : nei) {
                int v = knode.id;
                if (visited2[v] && dis2[v] + knode.dis < dis2[u2]) {
                    dis2[u2] = dis2[v] + knode.dis;
                    need_expand = 0;
                    break;
                }
            }
            if (visited1[u2]) {
                if (min_length > dis1[u2] + dis2[u2])
                    min_length = dis1[u2] + dis2[u2];

            }
            if (need_expand) {
                for (KNode knode : nei) {
                    int v = knode.id;
//					if (node_order[v] < node_order[u2]
//							|| (node_order[v] == node_order[u2]
//									&& inlev_order[v] < inlev_order[u2])
//							|| (node_order[v] == node_order[u2]
//									&& inlev_order[v] == inlev_order[u2]
//									&& v < u2)) {
                    if (node_order[v] < node_order[u2]) {
                        continue;
                    }
                    if (visited2[v]) {
                        if (knode.dis < dis2[v] - dis2[u2]) {
                            dis2[v] = dis2[u2] + knode.dis;
                            q2->decrease_priority(v, dis2[v]);
                        }
                    } else {
                        dis2[v] = dis2[u2] + knode.dis;
                        q2->push(v, dis2[v]);
                        visited2[v] = 1;
                        visitedVector2.push_back(v);
                    }
                }
            }
        }

    }
    for (int i = 0; i < visitedVector1.size(); i++) {
        dis1[visitedVector1[i]] = INT_MAX;
        visited1[visitedVector1[i]] = 0;
    }
    for (int i = 0; i < visitedVector2.size(); i++) {
        dis2[visitedVector2[i]] = INT_MAX;
        visited2[visitedVector2[i]] = 0;
    }
    q1->clear();
    q2->clear();
    return min_length;
}

int hierSP_rev(int s, int t, int *node_order, unordered_map<int, vector<KNode> > &shortcuts, mem_struct &mems,
               int thresh_dist) {
    int min_length = INT_MAX;
    DijkstraQueue *q1 = mems.q;
    DijkstraQueue *q2 = mems.q1;
    long *dis1 = mems.dist;
    long *dis2 = mems.dist1;
    int *visited1 = mems.visited;
    int *visited2 = mems.visited1;
    q1->clear();
    q2->clear();
    vector<int> visitedVector1;
    vector<int> visitedVector2;
    int valid1 = 1;
    int valid2 = 1;
    dis1[s] = 0;
    dis2[t] = 0;
    q1->push(s, 0);
    q2->push(t, 0);
    visited1[s] = 1;
    visited2[t] = 1;
    visitedVector1.push_back(s);
    visitedVector2.push_back(t);
//	cout<<"there"<<endl;
    while ((!q1->empty() && valid1) || (!q2->empty() && valid2)) {
        if (!q1->empty() && valid1) {
            int u1 = q1->pop();
            if (dis1[u1] > min_length || dis1[u1] > thresh_dist) {
                valid1 = 0;
            }
            visited1[u1] = 2;

            int need_expand = 1;
            vector<KNode> nei = shortcuts_inc_arr[u1];
            for (KNode knode : nei) {
                int v = knode.id;
                if (visited1[v] && dis1[v] + knode.dis < dis1[u1]) {
                    dis1[u1] = dis1[v] + knode.dis;
                    need_expand = 0;
                    break;
                }
            }
            if (visited2[u1]) {
                if (min_length > dis1[u1] + dis2[u1])
                    min_length = dis1[u1] + dis2[u1];

            }

            if (need_expand) {
                for (KNode knode : nei) {
                    int v = knode.id;
//					if (node_order[v] < node_order[u1]
//							|| (node_order[v] == node_order[u1]
//									&& inlev_order[v] < inlev_order[u1])
//							|| (node_order[v] == node_order[u1]
//									&& inlev_order[v] == inlev_order[u1]
//									&& v < u1)) {
                    if (node_order[v] < node_order[u1]) {
                        continue;
                    }
                    if (visited1[v]) {
                        if (knode.dis < dis1[v] - dis1[u1]) {
                            dis1[v] = dis1[u1] + knode.dis;
                            q1->decrease_priority(v, dis1[v]);
                        }
                    } else {
                        dis1[v] = dis1[u1] + knode.dis;
                        q1->push(v, dis1[v]);
                        visited1[v] = 1;
                        visitedVector1.push_back(v);

                    }

                }
            }
        }

        if (!q2->empty() && valid2) {
            int u2 = q2->pop();
            if (dis2[u2] > min_length || dis2[u2] > thresh_dist) {
                valid2 = 0;
            }
            visited2[u2] = 2;

            int need_expand = 1;
            vector<KNode> nei = shortcuts_arr[u2];
            for (KNode knode : nei) {
                int v = knode.id;
                if (visited2[v] && dis2[v] + knode.dis < dis2[u2]) {
                    dis2[u2] = dis2[v] + knode.dis;
                    need_expand = 0;
                    break;
                }
            }
            if (visited1[u2]) {
                if (min_length > dis1[u2] + dis2[u2])
                    min_length = dis1[u2] + dis2[u2];

            }
            if (need_expand) {
                for (KNode knode : nei) {
                    int v = knode.id;
//					if (node_order[v] < node_order[u2]
//							|| (node_order[v] == node_order[u2]
//									&& inlev_order[v] < inlev_order[u2])
//							|| (node_order[v] == node_order[u2]
//									&& inlev_order[v] == inlev_order[u2]
//									&& v < u2)) {
                    if (node_order[v] < node_order[u2]) {
                        continue;
                    }
                    if (visited2[v]) {
                        if (knode.dis < dis2[v] - dis2[u2]) {
                            dis2[v] = dis2[u2] + knode.dis;
                            q2->decrease_priority(v, dis2[v]);
                        }
                    } else {
                        dis2[v] = dis2[u2] + knode.dis;
                        q2->push(v, dis2[v]);
                        visited2[v] = 1;
                        visitedVector2.push_back(v);
                    }
                }
            }
        }

    }
    for (int i = 0; i < visitedVector1.size(); i++) {
        dis1[visitedVector1[i]] = INT_MAX;
        visited1[visitedVector1[i]] = 0;
    }
    for (int i = 0; i < visitedVector2.size(); i++) {
        dis2[visitedVector2[i]] = INT_MAX;
        visited2[visitedVector2[i]] = 0;
    }
    q1->clear();
    q2->clear();
    return min_length;
}

int test_shortest_path_bidirection(int s, int t, int *node_order, unordered_map<int, vector<KNode> > &shortcuts,
                                   mem_struct &mems,
                                   int thresh_dist) {
    int min_length = INT_MAX;
    DijkstraQueue *q1 = mems.q;
    DijkstraQueue *q2 = mems.q1;
    long *dis1 = mems.dist;
    long *dis2 = mems.dist1;
    int *visited1 = mems.visited;
    int *visited2 = mems.visited1;
    q1->clear();
    q2->clear();
    vector<int> visitedVector1;
    vector<int> visitedVector2;
    int *pre1 = new int[test_n + 1];
    int *pre2 = new int[test_n + 1];
    int valid1 = 1;
    int valid2 = 1;
    dis1[s] = 0;
    dis2[t] = 0;
    pre1[s] = 0;
    pre2[t] = 0;
    q1->push(s, 0);
    q2->push(t, 0);
    visited1[s] = 1;
    visited2[t] = 1;
    visitedVector1.push_back(s);
    visitedVector2.push_back(t);
//	cout<<"there"<<endl;
    int mid_node = 0;
    while ((!q1->empty() && valid1) || (!q2->empty() && valid2)) {
        if (!q1->empty() && valid1) {
            int u1 = q1->pop();
            if (dis1[u1] > min_length || dis1[u1] > thresh_dist) {
                valid1 = 0;

            }
            visited1[u1] = 2;

            int need_expand = 1;
            vector<KNode> nei = shortcuts_arr[u1];
//			for (KNode knode : nei) {
//				int v = knode.id;
//				if (visited1[v] && dis1[v] + knode.dis < dis1[u1]) {
//					dis1[u1] = dis1[v] + knode.dis;
//					need_expand = 0;
//					break;
//				}
//			}
            if (visited2[u1]) {
                if (min_length > dis1[u1] + dis2[u1]) {
                    min_length = dis1[u1] + dis2[u1];
                    mid_node = u1;
                }

            }

            if (need_expand) {
                for (KNode knode : nei) {
                    int v = knode.id;
//					if (node_order[v] < node_order[u1]
//							|| (node_order[v] == node_order[u1]
//									&& inlev_order[v] < inlev_order[u1])
//							|| (node_order[v] == node_order[u1]
//									&& inlev_order[v] == inlev_order[u1]
//									&& v < u1)) {
                    if (node_order[v] < node_order[u1]) {
                        continue;
                    }
                    if (visited1[v]) {
                        if (knode.dis < dis1[v] - dis1[u1]) {
                            dis1[v] = dis1[u1] + knode.dis;
                            pre1[v] = u1;
                            q1->decrease_priority(v, dis1[v]);
                        }
                    } else {
                        dis1[v] = dis1[u1] + knode.dis;
                        q1->push(v, dis1[v]);
                        pre1[v] = u1;
                        visited1[v] = 1;
                        visitedVector1.push_back(v);

                    }

                }
            }
        }

        if (!q2->empty() && valid2) {
            int u2 = q2->pop();
            if (dis2[u2] > min_length || dis2[u2] > thresh_dist) {
                valid2 = 0;
//				int tmp = u2;
//				while (tmp) {
//					cout << tmp << " ";
//					tmp = pre1[tmp];
//				}
//				cout << endl << endl;
            }
            visited2[u2] = 2;

            int need_expand = 1;
            vector<KNode> nei = shortcuts_arr[u2];
//			for (KNode knode : nei) {
//				int v = knode.id;
//				if (visited2[v] && dis2[v] + knode.dis < dis2[u2]) {
//					dis2[u2] = dis2[v] + knode.dis;
//					need_expand = 0;
//					break;
//				}
//			}
            if (visited1[u2]) {
                if (min_length > dis1[u2] + dis2[u2]) {
                    min_length = dis1[u2] + dis2[u2];
                    mid_node = u2;
                }

            }
            if (need_expand) {
                for (KNode knode : nei) {
                    int v = knode.id;
//					if (node_order[v] < node_order[u2]
//							|| (node_order[v] == node_order[u2]
//									&& inlev_order[v] < inlev_order[u2])
//							|| (node_order[v] == node_order[u2]
//									&& inlev_order[v] == inlev_order[u2]
//									&& v < u2)) {
                    if (node_order[v] < node_order[u2]) {
                        continue;
                    }
                    if (visited2[v]) {
                        if (knode.dis < dis2[v] - dis2[u2]) {
                            dis2[v] = dis2[u2] + knode.dis;
                            pre2[v] = u2;
                            q2->decrease_priority(v, dis2[v]);
                        }
                    } else {
                        dis2[v] = dis2[u2] + knode.dis;
                        q2->push(v, dis2[v]);
                        pre2[v] = u2;
                        visited2[v] = 1;
                        visitedVector2.push_back(v);
                    }
                }
            }
        }

    }
    for (int i = 0; i < visitedVector1.size(); i++) {
        dis1[visitedVector1[i]] = INT_MAX;
        visited1[visitedVector1[i]] = 0;
    }
    for (int i = 0; i < visitedVector2.size(); i++) {
        dis2[visitedVector2[i]] = INT_MAX;
        visited2[visitedVector2[i]] = 0;
    }
    int tmp = mid_node;
    vector<int> path;
//	cout << "paths in bilevel" << endl;
//	while (tmp) {
//		path.insert(path.begin(), tmp);
////						cout << tmp << " ";
//		tmp = pre1[tmp];
//	}
//	tmp = mid_node;
//	while (tmp) {
//		path.push_back(tmp);
//		//						cout << tmp << " ";
//		tmp = pre2[tmp];
//	}
//	for (int node : path) {
//		cout << node << " " << node_order[node] << endl;
//	}
//	cout << endl;
    q1->clear();
    q2->clear();
    delete[] pre1;
    delete[] pre2;
    return min_length;
}

int hierSP_level(int s, int t, vector<int> &path, int level) {
    int *dis = new int[test_n + 1];
    int *pre = new int[test_n + 1];
    DijkstraQueue *q = new DijkstraQueue(test_n + 1);
    int *visited = new int[test_n + 1];
    memset(visited, 0, sizeof(int) * (test_n + 1));
    vector<int> visitedVector;
    dis[s] = 0;
    visitedVector.push_back(s);
    visited[s] = 1;
    pre[s] = 0;
    q->clear();
    q->push(s, 0);
    while (!q->empty()) {
        int u = q->pop();
        if (u == t) {
            int temp = t;
            do {
                path.insert(path.begin(), temp);
                temp = pre[temp];
            } while (temp != 0);
            int r = dis[u];

            return r;
        }
        for (KNode knode : level_scs_arr[level][u]) {
            int v = knode.id;

            if (visited[v]) {
                if (knode.dis < dis[v] - dis[u]) {
                    dis[v] = dis[u] + knode.dis;
                    pre[v] = u;
                    q->decrease_priority(v, dis[v]);
                }
            } else {
                dis[v] = dis[u] + knode.dis;
                q->push(v, dis[v]);
                pre[v] = u;
                visited[v] = 1;
                visitedVector.push_back(v);
            }

        }
    }
    for (int node : visitedVector) {
        visited[node] = 0;
    }
    delete q;
    delete[] dis;
    delete[] pre;
    delete[] visited;
    return 0;
}




void get_distance(int s, int * dis) {
    for (int i = 1; i < test_n + 1; ++i)
        dis[i] = INT_MAX;
    dis[s] = 0;
    DijkstraQueue *q = new DijkstraQueue(test_n, s);
    while (!q->empty()) {
        int u = q->pop();
        for (int i = vStart[u]; i != 0; i = eNext[i]) {
            int v = eNode[i];
            if (eDis[i] < dis[v] - dis[u]) {
                dis[v] = dis[u] + eDis[i];
                q->decrease_priority(v, dis[v]);
            }
        }
    }
    delete q;
}



int get_distance_st(int s, int t) {
    int* dis = new int[test_n + 1];
    int* pre = new int[test_n + 1];
    for (int i = 1; i < test_n + 1; ++i)
        dis[i] = INT_MAX;
    dis[s] = 0;
    pre[s] = 0;
    DijkstraQueue *q = new DijkstraQueue(test_n, s);
    //	int sx = cell_loc[s] / cell_num, sy = cell_loc[s] % cell_num;
    while (!q->empty()) {
        int u = q->pop();
        if (u == t) {
            int r = dis[u];
            delete q;
            delete[] dis;
            delete[] pre;
            return r;
        }
        for (int i = vStart[u]; i != 0; i = eNext[i]) {
            int v = eNode[i];
            if (eDis[i] < dis[v] - dis[u]) {
                dis[v] = dis[u] + eDis[i];
                pre[v] = u;
                q->decrease_priority(v, dis[v]);
            }
        }
    }
    delete q;
    delete[] dis;
    delete[] pre;
    return 0;
}

int get_path_st(int s, int t, vector<int>& path) {
    int* dis = new int[test_n + 1];
    int* pre = new int[test_n + 1];
    DijkstraQueue* q = new DijkstraQueue(test_n + 1);
    for (int i = 1; i < test_n + 1; ++i)
        dis[i] = INT_MAX;
    unordered_map<int, int> tested_nodes;
    tested_nodes[s] = 1;
    dis[s] = 0;
    pre[s] = 0;
    q->clear();
    q->init(s);
    while (!q->empty()) {
        int u = q->pop();
        if (u == t) {
            int temp = t;
            do {
                path.insert(path.begin(), temp);
                temp = pre[temp];
            } while (temp != 0);
            int r = dis[u];

            return r;
        }
        for (int i = vStart[u]; i != 0; i = eNext[i]) {
            int v = eNode[i];

            if (eDis[i] < dis[v] - dis[u]) {
                dis[v] = dis[u] + eDis[i];
                pre[v] = u;
                q->decrease_priority(v, dis[v]);
            }

        }
    }
    delete q;
    delete[] dis;
    delete[] pre;
    return 0;
}
#endif //SOB_TESTSHORTESTPATH_H
