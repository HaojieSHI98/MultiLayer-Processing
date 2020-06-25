/*
 * HierIndex.h
 *
 *  Created on: 2016年9月5日
 *      Author: lostrong
 */

#ifndef SRC_HIERINDEX_H_
#define SRC_HIERINDEX_H_

#include<vector>
#include<limits.h>

#include "../util/MemoryAllocation.h"
#include "../dijk/NaiveKNN.h"
#include "../util/GeneralUtility.h"
#include "../testing/TestShortestPath.h"
#include "../io/Input.h"
#include "ScobConfiguration.h"
#include "../RoadKNNUtilities.h"

using namespace std;
// num_nodes is the maximal number of nodes allowed in a cell

int build_cells() {
    x_min = INT_MAX, y_min = INT_MAX, x_max = INT_MIN, y_max = INT_MIN, x_len = 0, y_len = 0;
    Point *a = new Point[test_n + 1];
    vector<Point> *b = new vector<Point>[cell_num];
    cout << test_n << endl;
    for (int i = 0; i < test_n; ++i) {
        x_min = min(x_min, x[i]);
        x_max = max(x_max, x[i]);
        y_min = min(y_min, y[i]);
        y_max = max(y_max, y[i]);
    }
    int gap = (y_max - y_min) - (x_max - x_min);
    cout << "gap:" << gap << endl;
    cout << "x_min:" << x_min << "x_max:" << x_max << "(y_max+y_min)/2:" << (y_max + y_min) / 2 << "(y_max + y_min) / 2"
         << (y_max + y_min) / 2 << endl;
    largest_distance = get_distance_lat_lng(x_min, (y_max + y_min) / 2, x_max, (y_max + y_min) / 2, 1000000.0);
    cout << "largest_distance:" << largest_distance << endl;
    if (gap > 0) {
        //extend x_len
        x_min -= gap / 2;
        x_max += gap / 2;
    } else {
        gap = -gap;
        y_min -= gap / 2;
        y_max += gap / 2;
    }
    x_len = (x_max - x_min) / cell_num + 1;
    y_len = (y_max - y_min) / cell_num + 1;

    for (int i = 1; i < test_n + 1; ++i) {
        a[i].x = x[i - 1], a[i].y = y[i - 1];
        a[i].flag = i;
    }
    sort(a + 1, a + test_n + 1, Point::compare_x);
    for (int i = 1; i < test_n + 1; ++i) {
        int j = (a[i].x - x_min) / x_len;
        b[j].push_back(a[i]);
    }
    delete a;
//#pragma omp parallel for
    for (int i = 0; i < cell_num; ++i) {

        sort(b[i].begin(), b[i].end(), Point::compare_y);
        for (int j = 0; j < b[i].size(); ++j) {
            int k = (b[i][j].y - y_min) / y_len;
            cells[i][k].push_back(b[i][j].flag);
            cell_loc[b[i][j].flag] = i * cell_num + k;
            cell_loc_x[b[i][j].flag] = i;
            cell_loc_y[b[i][j].flag] = k;
            if (k >= cell_num) {
                cout << "k>=cell_num" << endl;
                while (true);
            }
        }
    }
    const int xc[25] = {-2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2};
    const int yc[25] = {-2, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};

#pragma omp parallel for
    for (int i = 0; i < cell_num; ++i) {

        for (int j = 0; j < cell_num; ++j) {
            for (int k = 0; k < 25; ++k) {
                int ux = i + xc[k], uy = j + yc[k];
                if (ux >= 0 && ux < cell_num && uy >= 0 && uy < cell_num) {
                    size_t len = circle_cells[i][j].size();
                    for (int l = 0; l < cells[ux][uy].size(); ++l) {
                        circle_cell_loc[(xc[k] + 2) * 5 + yc[k] + 2][cells[ux][uy][l]] = len + l;
                        int loc = cell_loc[cells[ux][uy][l]];
                        if (loc / cell_num != ux || loc % cell_num != uy) {
                            cout << "asign false" << endl;
                            cout << ux << " " << uy << " " << l << endl;
                            while (true);
                        }

                        int ts = get_location(i, j, cells[ux][uy][l]);
                        if (ts != len + l) {
                            cout << "inconsistent!" << endl;
                            while (true);
                        }

                    }

                    for (auto it = cells[ux][uy].begin(); it != cells[ux][uy].end(); it++)
                        circle_cells[i][j].insert(circle_cells[i][j].end(), *it);

                }
            }
        }
    }
    cell_loc[0] = -1;

    return 0;
}

int get_cellnum_by_number(int num_nodes) {
    int div_len = INT_MAX;
    build_cells();
    queue<CellNode> que;
    //init queue
    int total = 0;
    for (int i = 0; i < cell_num; ++i) {
        for (int j = 0; j < cell_num; ++j) {
            CellNode node(i * x_len + x_min, (i + 1) * x_len + x_min, j * y_len + y_min, (j + 1) * y_len + y_min,
                          cells[i][j]);
            que.push(node);
            total += cells[i][j].size();
        }
    }

    cout << "total nodes: " << total << endl;
    while (!que.empty()) {
        CellNode front_node = que.front();
        if (front_node.inside_nodes.size() > num_nodes) {
            int x_min = front_node.coor_x_min;
            int x_max = front_node.coor_x_max;
            int x_mid = (front_node.coor_x_min + front_node.coor_x_max) / 2;

            int y_min = front_node.coor_y_min;
            int y_max = front_node.coor_y_max;
            int y_mid = (front_node.coor_y_min + front_node.coor_y_max) / 2;

            vector<int> inside_nodes1, inside_nodes2, inside_nodes3, inside_nodes4;
            for (int i = 0; i < front_node.inside_nodes.size(); i++) {
                int node_id = front_node.inside_nodes[i];
                int coor_x = x[node_id - 1];
                int coor_y = y[node_id - 1];
                if (coor_x <= x_mid && coor_x >= x_min && coor_y <= y_mid && coor_y >= y_min)
                    inside_nodes1.push_back(node_id);
                if (coor_x <= x_mid && coor_x >= x_min && coor_y >= y_mid && coor_y <= y_max)
                    inside_nodes2.push_back(node_id);
                if (coor_x >= x_mid && coor_x <= x_max && coor_y <= y_mid && coor_y >= y_min)
                    inside_nodes3.push_back(node_id);
                if (coor_x >= x_mid && coor_x <= x_max && coor_y >= y_mid && coor_y <= y_max)
                    inside_nodes4.push_back(node_id);
            }

            CellNode node1(x_min, x_mid, y_min, y_mid, inside_nodes1);
            CellNode node2(x_min, x_mid, y_mid, y_max, inside_nodes2);
            CellNode node3(x_mid, x_max, y_min, y_mid, inside_nodes3);
            CellNode node4(x_mid, x_max, y_mid, y_max, inside_nodes4);
            que.pop();
            que.push(node1);
            que.push(node2);
            que.push(node3);
            que.push(node4);
            //		cout<<"expand"<<endl;
        } else {

            if (div_len > front_node.coor_x_max - front_node.coor_x_min)
                div_len = front_node.coor_x_max - front_node.coor_x_min;
            que.pop();

        }

    }

    return (x_max - x_min) / (div_len + 1);
}

void cal_node_order_level(int cell_num, int level) {
    cout << "start init cells..." << endl;
    init_cells(cell_num);
    build_cells();
    cout << "start getting landmarks..." << endl;
    get_regional_landmarks(0, 1, level - 1);
    for (int i = 0; i < test_n; ++i)
        if (regional_landmarks_all[i])
            node_order[i] = level;
    memset(regional_landmarks_all, 0, sizeof(int) * test_n);
    dispose_cells(cell_num);
}

int cal_node_order(int num_nodes) {
    cell_num = 10;
    init_cells(cell_num);
    int temp_num = get_cellnum_by_number(num_nodes);
    dispose_cells(cell_num);
    cout << "temp_num: " << temp_num << endl;
    cell_num = 1;
    int level_ratio = 2;
    while (temp_num > 1) {
        temp_num /= level_ratio;
        cell_num *= level_ratio;
    }
    int level = 1;
    while (cell_num > 5 && level <= max_level) {
        cout << "level: " << level << endl;
        cout << "cell_num: " << cell_num << endl;
        cell_nums[level] = cell_num;
        cal_node_order_level(cell_num, level);
        level++;
        cell_num /= level_ratio;
    }
    cout << "the number of levels: " << level << endl;
    return level;
}


int cal_node_order_coverage(int num_nodes) {
    FILE* file = fopen((input_parameters.input_data_dir+"1coverage_node_order.txt").c_str(), "r");
    int nodeid, level;
    while(fscanf(file, "%d %d\n",&nodeid, &level)==2){
        node_order[nodeid]=level;
    }
    fclose(file);
    return 7;
}


int cal_node_order_between(int num_nodes) {
    FILE* file = fopen((input_parameters.input_data_dir+"1between_node_order.txt").c_str(), "r");
    int nodeid, level;
    while(fscanf(file, "%d %d\n",&nodeid, &level)==2){
        node_order[nodeid]=level;
    }
    fclose(file);
    return 7;
}


//depreciated
unordered_map<int, vector<KNode> > add_shortcuts() {
    string file_name1 = data_space;
    file_name1 = file_name1.append("hiershortcuts.txt");
    FILE *shortcut_file = fopen(file_name1.c_str(), "w");
    unordered_map<int, vector<KNode> > shortcuts;
    int *visited = new int[test_n + 1];
    memset(visited, 0, sizeof(int) * (test_n + 1));
    int *dis = new int[test_n + 1];
    int *marked_first = new int[test_n + 1];
    int *pre = new int[test_n + 1];
    memset(marked_first, 0, sizeof(int) * (test_n + 1));
    DijkstraQueue *q = new DijkstraQueue(test_n + 1);
    for (int i = 1; i <= test_n; ++i) {
        int query = i;
        vector<int> visitedlist;
        //store shortcuts
        vector<KNode> result;
        dis[query] = 0;
        pre[query] = 0;
        visited[query] = 1;
        visitedlist.push_back(query);
        q->push(query, 0);
        int rank = node_order[query];
//		int inlev_rank = inlev_order[query];
        int unmarked_cnt = 1;
        while (!q->empty()) {
            int u = q->pop();
            if (marked_first[u] == 0) {
                unmarked_cnt--;
            }
//			if ((node_order[u] > rank
//					|| (node_order[u] == rank && inlev_order[u] > inlev_rank))
//					&& marked_first[u] == 0) {
            if (node_order[u] >= rank && marked_first[u] == 0) {
                if (u != query) {
                    result.push_back(KNode(u, dis[u]));

                    fprintf(shortcut_file, "%d %d %d\n", query, u, dis[u]);
                    marked_first[u] = 1;
                }

            }
            for (int i = vStart[u]; i != 0; i = eNext[i]) {
                int v = eNode[i];
                if (visited[v]) {
                    if (eDis[i] < dis[v] - dis[u]) {
                        dis[v] = dis[u] + eDis[i];
                        pre[v] = u;
                        if (marked_first[v] == 0)
                            unmarked_cnt--;
                        if (marked_first[u] == 0)
                            unmarked_cnt++;
                        marked_first[v] = marked_first[u];
                        q->decrease_priority(v, dis[v]);
                    }
                } else {
                    dis[v] = dis[u] + eDis[i];
                    q->push(v, dis[v]);
                    pre[v] = u;
                    marked_first[v] = marked_first[u];
                    if (marked_first[v] == 0)
                        unmarked_cnt++;
                    visited[v] = 1;
                    visitedlist.push_back(v);
                }

            }
            if (unmarked_cnt == 0) {
                break;
            }
        }

        for (int item : visitedlist) {
            visited[item] = 0;
            marked_first[item] = 0;
        }

        q->clear();
        //	cout << "delete time: " << clock() - start << endl;
        shortcuts[query] = result;
    }
    fclose(shortcut_file);
    delete[] visited;
    delete[] dis;
    delete[] marked_first;
    delete[] pre;
    delete q;
    return shortcuts;
}


void build_scob_one_level(int level, FILE *shortcut_file, FILE *shortcut_inc_file, FILE *shortcut_file_hori,
                          FILE **level_scs,
                          vector<KNode> *shortcuts_tmp_arr_1, vector<KNode> *shortcuts_tmp_arr_2) {
    int *visited = new int[test_n + 1];
    memset(visited, 0, sizeof(int) * (test_n + 1));
    int *dis = new int[test_n + 1];
    int *marked_first = new int[test_n + 1];
    vector<int> *pre = new vector<int>[test_n + 1];
    memset(marked_first, 0, sizeof(int) * (test_n + 1));
    DijkstraQueue *q = new DijkstraQueue(test_n + 1);
    for (int i = 1; i <= test_n; ++i) {
//        cout<<i<<endl;
        int query = i;
        int rank = node_order[query];
        if (rank < level - 1)
            continue;
        vector<int> visitedlist;
        //store shortcuts
        vector<KNode> result;
        dis[query] = 0;
        pre[query].push_back(0);
        visited[query] = 1;
        visitedlist.push_back(query);
        q->push(query, 0);

//		int inlev_rank = inlev_order[query];
        int unmarked_cnt = 1;
        int num_nodes_visited=0;
        while (!q->empty()) {
            num_nodes_visited++;
            int u = q->pop();
            if (marked_first[u] == 0) {
                unmarked_cnt--;
            }
            if (node_order[u] >= level) {
                if (u != query && rank == level) {
                    fprintf(shortcut_file_hori, "%d %d %d\n", query, u, dis[u]);

                }
            }
            if (node_order[u] >= level && marked_first[u] == 0) {
                if (u != query && rank == level) {
                    fprintf(shortcut_file, "%d %d %d\n", query, u, dis[u]);

                }
                if (u != query && rank == level - 1) {
                    fprintf(shortcut_inc_file, "%d %d %d\n", query, u, dis[u]);
                }

                if (u != query && rank >= level) {
                    result.push_back(KNode(u, dis[u]));

                }
                if (u != query) {
                    marked_first[u] = 1;
                }

            }
            for (KNode knode : shortcuts_tmp_arr_1[u]) {
                int v = knode.id;
                int disv = knode.dis;
                if (visited[v]) {
                    if (disv < dis[v] - dis[u]) {
                        dis[v] = dis[u] + disv;
                        pre[v].clear();
                        pre[v].push_back(u);
                        if (marked_first[v] == 0)
                            unmarked_cnt--;
                        if (marked_first[u] == 0)
                            unmarked_cnt++;
                        marked_first[v] = marked_first[u];
                        q->decrease_priority(v, dis[v]);
                    } else if (disv == dis[v] - dis[u]) {
                        pre[v].push_back(u);
                        if (marked_first[v] == 0)
                            unmarked_cnt--;
                        marked_first[v] = 1;
                        for (int pre_node : pre[v]) {
                            if (marked_first[pre_node] == 0) {
                                marked_first[v] = 0;
                                unmarked_cnt++;
                                break;
                            }
                        }
                    }
                } else {
                    dis[v] = dis[u] + disv;
                    q->push(v, dis[v]);
                    pre[v].push_back(u);
                    marked_first[v] = marked_first[u];
                    if (marked_first[v] == 0)
                        unmarked_cnt++;
                    visited[v] = 1;
                    visitedlist.push_back(v);
                }

            }

            if (unmarked_cnt == 0) {
                break;
            }

        }
//        cout<<"visited nodes: "<<num_nodes_visited<<endl;
        for (int item : visitedlist) {
            visited[item] = 0;
            marked_first[item] = 0;
        }

        q->clear();

        for (KNode knode : result) {
            shortcuts_tmp_arr_2[query].push_back(KNode(knode.id, knode.dis));
            fprintf(level_scs[level], "%d %d %d\n", query, knode.id, knode.dis);
        }
    }
    for (int i = 0; i <= test_n; i++) {
        shortcuts_tmp_arr_1[i].clear();
        pre[i].clear();
    }
    delete[] visited;
    delete[] dis;
    delete[] marked_first;
    delete[] pre;
    delete q;

}

void build_scob_level_by_level(int highest_level) {
    cout<<"highest_level: "<<highest_level<<endl;
    cout << "start adding shortcuts level by level..." << endl;
    vector<KNode> *shortcuts_tmp_arr_1 = new vector<KNode>[MAGIC_NUM];
    vector<KNode> *shortcuts_tmp_arr_2 = new vector<KNode>[MAGIC_NUM];

    cout << "starting creating files..." << endl;
    string file_name1 = data_space;
    if(COVERAGE_INDEX)
        file_name1=file_name1.append("coverage_hiershortcuts.txt");
    else if(BETWEEN_INDEX)
        file_name1=file_name1.append("between_hiershortcuts.txt");
    else if(ADA_COVERAGE_INDEX)
        file_name1=file_name1.append("ada_coverage_hiershortcuts.txt");
    else if(ADA_BETWEEN_INDEX)
        file_name1=file_name1.append("ada_between_hiershortcuts.txt");
    else
        file_name1 = file_name1.append("hiershortcuts.txt");
    file_name1 = input_parameters.input_data_dir+file_name1;
    cout<<"shortcut file name: "<<file_name1<<endl;
    FILE *shortcut_file = fopen(file_name1.c_str(), "w");

    string file_name2 = data_space;
    if(COVERAGE_INDEX)
        file_name2=file_name2.append("hiershortcuts_strict_inc.txt");
    else if(BETWEEN_INDEX)
        file_name2=file_name2.append("hiershortcuts_strict_inc.txt");
    else if(ADA_COVERAGE_INDEX)
        file_name2=file_name2.append("ada_coverage_hiershortcuts_strict_inc.txt");
    else if(ADA_BETWEEN_INDEX)
        file_name2=file_name2.append("ada_between_hiershortcuts_strict_inc.txt");
    else
        file_name2 = file_name2.append("hiershortcuts_strict_inc.txt");
    file_name2 = input_parameters.input_data_dir+file_name2;
    FILE *shortcut_file_inc = fopen(file_name2.c_str(), "w");

    string file_name3 = data_space;
    if(COVERAGE_INDEX)
        file_name3=file_name3.append("hiershortcuts_strict_hori.txt");
    else if(BETWEEN_INDEX)
        file_name3=file_name3.append("hiershortcuts_strict_hori.txt");
    else if(ADA_COVERAGE_INDEX)
        file_name3=file_name3.append("ada_coverage_hiershortcuts_strict_hori.txt");
    else if(ADA_BETWEEN_INDEX)
        file_name3=file_name3.append("ada_between_hiershortcuts_strict_hori.txt");
    else
        file_name3 = file_name3.append("hiershortcuts_strict_hori.txt");
    file_name3 = input_parameters.input_data_dir+file_name3;
    FILE *shortcut_file_hori = fopen(file_name3.c_str(), "w");

    FILE **level_scs = new FILE *[highest_level];
    for (int i = 0; i <= highest_level; i++) {
        string file_name = data_space;
        if(COVERAGE_INDEX)
            file_name=file_name.append("coverage_level_sc");
        else if(BETWEEN_INDEX)
            file_name=file_name.append("between_level_sc");
        else if(ADA_COVERAGE_INDEX)
            file_name2=file_name2.append("ada_coverage_level_sc.txt");
        else if(ADA_BETWEEN_INDEX)
            file_name2=file_name2.append("ada_between_level_sc.txt");
        else
            file_name = file_name.append("level_sc");
        stringstream ss;
        ss << i;
        string str = ss.str();
        file_name = file_name.append(str);
        file_name = input_parameters.input_data_dir+file_name;
        cout<<"level file name: "<<file_name<<endl;
        level_scs[i] = fopen(file_name.c_str(), "w");
    }

    // add level 0 shortcuts, e.g., original edges
    cout << "starting level 0 shortcuts..." << endl;
    for (int u = 0; u <= test_n; u++) {
        for (int i = vStart[u]; i != 0; i = eNext[i]) {
            shortcuts_tmp_arr_1[u].push_back(KNode(eNode[i], eDis[i]));
            if (node_order[u] == 0) {
                fprintf(shortcut_file, "%d %d %d\n", u, eNode[i], eDis[i]);

            }
            fprintf(level_scs[0], "%d %d %d\n", u, eNode[i], eDis[i]);
        }
    }

    for (int level = 1; level <= highest_level; level++) {
        cout << "starting level " << level << "..." << endl;
        if (level % 2 == 1) {
            build_scob_one_level(level, shortcut_file, shortcut_file_inc, shortcut_file_hori, level_scs,
                                 shortcuts_tmp_arr_1,
                                 shortcuts_tmp_arr_2);
        } else {
            build_scob_one_level(level, shortcut_file, shortcut_file_inc, shortcut_file_hori, level_scs,
                                 shortcuts_tmp_arr_2,
                                 shortcuts_tmp_arr_1);

        }
        cout<<"ends level "<<level<<endl;

    }
    fclose(shortcut_file);
    cout<<"shortcut_file closed"<<endl;
    fclose(shortcut_file_inc);
    cout<<"shortcut_file_inc closed"<<endl;
    fclose(shortcut_file_hori);
    cout<<"shortcut_file_hori closed"<<endl;

    for (int i = 0; i <= highest_level; i++) {
        fclose(level_scs[i]);
        cout<<"level_scs "<<i<<" closed"<<endl;
    }
    cout<<"files closed"<<endl;
    delete[] shortcuts_tmp_arr_1;
    delete[] shortcuts_tmp_arr_2;
    delete[] level_scs;
}

//deprecated
unordered_map<int, vector<KNode> > add_shortcuts_strict_inc() {
    string file_name1 = data_space;
    file_name1 = file_name1.append("hiershortcuts_strict_inc.txt");
    FILE *shortcut_file = fopen(file_name1.c_str(), "w");
    unordered_map<int, vector<KNode> > shortcuts;
    int *visited = new int[test_n + 1];
    memset(visited, 0, sizeof(int) * (test_n + 1));
    int *dis = new int[test_n + 1];
    int *marked_first = new int[test_n + 1];
    int *pre = new int[test_n + 1];
    memset(marked_first, 0, sizeof(int) * (test_n + 1));
    DijkstraQueue *q = new DijkstraQueue(test_n + 1);
    for (int i = 1; i <= test_n; ++i) {
        int query = i;
        vector<int> visitedlist;
        //store shortcuts
        vector<KNode> result;
        dis[query] = 0;
        pre[query] = 0;
        visited[query] = 1;
        visitedlist.push_back(query);
        q->push(query, 0);
        int rank = node_order[query];
        int inlev_rank = inlev_order[query];
        int unmarked_cnt = 1;
        while (!q->empty()) {
            int u = q->pop();
            if (marked_first[u] == 0) {
                unmarked_cnt--;
            }
            if (node_order[u] > rank && marked_first[u] == 0) {
                if (u != query) {
                    result.push_back(KNode(u, dis[u]));

                    fprintf(shortcut_file, "%d %d %d\n", query, u, dis[u]);
                    marked_first[u] = 1;
                }

            }
            for (int i = vStart[u]; i != 0; i = eNext[i]) {
                int v = eNode[i];
                if (visited[v]) {
                    if (eDis[i] < dis[v] - dis[u]) {
                        dis[v] = dis[u] + eDis[i];
                        pre[v] = u;
                        if (marked_first[v] == 0)
                            unmarked_cnt--;
                        if (marked_first[u] == 0)
                            unmarked_cnt++;
                        marked_first[v] = marked_first[u];
                        q->decrease_priority(v, dis[v]);
                    }
                } else {
                    dis[v] = dis[u] + eDis[i];
                    q->push(v, dis[v]);
                    pre[v] = u;
                    marked_first[v] = marked_first[u];
                    if (marked_first[v] == 0)
                        unmarked_cnt++;
                    visited[v] = 1;
                    visitedlist.push_back(v);
                }

            }
            if (unmarked_cnt == 0) {
                break;
            }
        }

        for (int item : visitedlist) {
            visited[item] = 0;
            marked_first[item] = 0;
        }

        q->clear();
        //	cout << "delete time: " << clock() - start << endl;
        shortcuts[query] = result;
    }
    fclose(shortcut_file);
    delete[] visited;
    delete[] dis;
    delete[] marked_first;
    delete[] pre;
    delete q;
    return shortcuts;
}

void build_scob() {
    cout<<"start building index..."<<endl;
    read_road_network();
    node_order = new int[test_n + 1];
    memset(node_order, 0, sizeof(int) * (test_n + 1));
    inlev_order = new int[test_n + 1];
    for (int u = 0; u < test_n + 1; u++) {
        inlev_order[u] = 0;
    }
//	cout << "start calculating node order..." << endl;
//    int highest_level = cal_node_order(leaf_num);
    int highest_level;
    if(COVERAGE_INDEX)
        highest_level= cal_node_order_coverage(leaf_num);
    else if(BETWEEN_INDEX)
        highest_level= cal_node_order_between(leaf_num);
    else highest_level = cal_node_order(leaf_num);

//	cout << "start adding shortcuts..." << endl;
    string file_name1 = data_space;
    if(COVERAGE_INDEX)
        file_name1=file_name1.append("coverage_node_order.txt");
    else if(BETWEEN_INDEX)
        file_name1=file_name1.append("between_node_order.txt");
    else if(ADA_COVERAGE_INDEX)
        file_name1=file_name1.append("ada_coverage_node_order.txt");
    else if(ADA_BETWEEN_INDEX)
        file_name1=file_name1.append("ada_between_node_order.txt");
    else
        file_name1 = file_name1.append("node_order.txt");
    file_name1 = input_parameters.input_data_dir+file_name1;
    FILE *nodeorder_file;
    nodeorder_file = fopen(file_name1.c_str(), "w");
    for (int i = 0; i < test_n + 1; ++i)
        fprintf(nodeorder_file, "%d %d\n", i, node_order[i]);
    fclose(nodeorder_file);

    nodeorder_file = fopen(file_name1.c_str(), "r");
    test_valid_file(nodeorder_file, file_name1, "build_scob");

    int node_id, rank;
    while (fscanf(nodeorder_file, "%d %d", &node_id, &rank) == 2) {
        node_order[node_id] = rank;
    }

    cout<<"finish ranking node orders..."<<endl;
//	add_shortcuts();
    string file_name2 = data_space;
    if(COVERAGE_INDEX)
        file_name2=file_name2.append("coverage_inlev_node_order.txt");
    else if(BETWEEN_INDEX)
        file_name2=file_name2.append("between_inlev_node_order.txt");
    else if(ADA_COVERAGE_INDEX)
        file_name2=file_name2.append("ada_coverage_inlev_node_order.txt");
    else if(ADA_BETWEEN_INDEX)
        file_name2=file_name2.append("ada_between_inlev_node_order.txt");
    else
        file_name2 = file_name2.append("inlev_node_order.txt");
    file_name2 = input_parameters.input_data_dir+file_name2;
    cout<<"inlev_node_order_name: "<<file_name2<<endl;
    FILE *inlev_file;
    inlev_file = fopen(file_name2.c_str(), "w");

    for (int i = 0; i < test_n; i++) {
        fprintf(inlev_file, "%d %d\n", i, inlev_order[i]);
    }
    fclose(inlev_file);
    inlev_file = fopen(file_name2.c_str(), "r");
    while (fscanf(inlev_file, "%d %d", &node_id, &rank) == 2) {
        inlev_order[node_id] = rank;
    }
    fclose(inlev_file);
    build_scob_level_by_level(highest_level);
//	add_shortcuts_strict_inc();
    delete[] node_order;
    delete[] inlev_order;
}


void read_scob_index() {
    if(has_read_scob_index==0 || memory_deleted==1) {
        string file_name1 = data_space;
        if(COVERAGE_INDEX)
            file_name1=file_name1.append("coverage_node_order.txt");
        else if(BETWEEN_INDEX)
            file_name1=file_name1.append("between_node_order.txt");
        else
            file_name1 = file_name1.append("node_order.txt");

        FILE *nodeorder_file = fopen((input_parameters.input_data_dir + file_name1).c_str(), "r");
        if (nodeorder_file == NULL) {
            cout << "nodeorder_file is NULL" << endl;
            cout<<"file name: "<<(input_parameters.input_data_dir + file_name1).c_str()<<endl;
            stop_here();
        }
        int node_id, rank;
        cout << "start reading node orders..." << endl;
        if (node_order == NULL)
            node_order = new int[test_n + 1];

        while (fscanf(nodeorder_file, "%d %d", &node_id, &rank) == 2) {
            node_order[node_id] = rank;
        }

        string file_name5 = data_space;
        if(COVERAGE_INDEX)
            file_name5=file_name5.append("coverage_inlev_node_order.txt");
        else if(BETWEEN_INDEX)
            file_name5=file_name5.append("between_inlev_node_order.txt");
        else
            file_name5 = file_name5.append("inlev_node_order.txt");
        FILE *inlev_nodeorder_file = fopen((input_parameters.input_data_dir + file_name5).c_str(), "r");
        test_valid_file(inlev_nodeorder_file, (input_parameters.input_data_dir + file_name5).c_str(),
                        "read_scob_index");

        cout << "start reading inlevel node orders..." << endl;
        if (inlev_order == NULL)
            inlev_order = new int[test_n + 1];

        memset(inlev_order, 0, sizeof(int) * (test_n + 1));
        cout << "start reading shortcuts..." << endl;
        string file_name2 = data_space;
        if(COVERAGE_INDEX)
            file_name2=file_name2.append("coverage_hiershortcuts.txt");
        else if(BETWEEN_INDEX)
            file_name2=file_name2.append("between_hiershortcuts.txt");
        else
            file_name2 = file_name2.append("hiershortcuts.txt");
        FILE *shortcut_file = fopen((input_parameters.input_data_dir + file_name2).c_str(), "r");
        cout<<"file name: "<<(input_parameters.input_data_dir + file_name2).c_str()<<endl;
        if (shortcut_file == NULL) {
            cout << "shortcut_file is NULL" << endl;
            stop_here();
        }
        int u, v, dis;
        int cnt_shortcuts = 0;
        while (fscanf(shortcut_file, "%d %d %d", &u, &v, &dis) == 3) {
            int flg = 0;
            for (KNode knode : shortcuts_arr[u]) {
                if (knode.id == v && knode.dis == dis) {
                    flg = 1;
                    break;
                }

            }
            if (flg == 0) {
                KNode *knode1 = new KNode(v, dis);
                shortcuts_arr[u].push_back(*knode1);
                cnt_shortcuts++;
                KNode *knode2 = new KNode(u, dis);
                rev_shortcuts_level_arr[v].push_back(*knode2);
            }

        }

        cout << "#shortcuts: " << cnt_shortcuts << endl;
        string file_name3 = data_space;
        //../ROADKNN-data/1
        if(COVERAGE_INDEX)
            file_name3=file_name3.append("coverage_hiershortcuts_strict_inc.txt");
        else if(BETWEEN_INDEX)
            file_name3=file_name3.append("between_hiershortcuts_strict_inc.txt");
        else
            file_name3 = file_name3.append("hiershortcuts_strict_inc.txt");
        FILE *shortcut_inc_file = fopen((input_parameters.input_data_dir + file_name3).c_str(), "r");
        if (shortcut_inc_file == NULL) {
            cout << "shortcut_inc_file is NULL" << endl;
            stop_here();
        }
        int cnt_shortcuts_inc = 0;
        while (fscanf(shortcut_inc_file, "%d %d %d", &u, &v, &dis) == 3) {
            int flg = 0;
            // test repeat shortcuts
            for (KNode knode : shortcuts_inc_arr[u]) {
                if (knode.id == v && knode.dis == dis) {
                    flg = 1;
                    break;
                }

            }
            if (flg == 0) {
                KNode *c = new KNode(v, dis);
                shortcuts_inc_arr[u].push_back(*c);
                cnt_shortcuts_inc++;
                KNode *knode4 = new KNode(u, dis);
                rev_shortcuts_arr[v].push_back(*knode4);

            }

        }

        cout << "#shortcuts_inc: " << cnt_shortcuts_inc << endl;

        fclose(nodeorder_file);
        fclose(inlev_nodeorder_file);
        fclose(shortcut_file);
        fclose(shortcut_inc_file);


        int lv_num = get_level_num_with_dataspace();
        FILE **level_scs = new FILE *[lv_num];

        for (int i = 0; i <= lv_num; i++) {
            string file_name = data_space;
            if(COVERAGE_INDEX)
                file_name=file_name.append("coverage_level_sc");
            else if(BETWEEN_INDEX)
                file_name=file_name.append("between_level_sc");
            else
                file_name = file_name.append("level_sc");
            stringstream ss;
            ss << i;
            string str = ss.str();
            file_name = file_name.append(str);
            level_scs[i] = fopen((input_parameters.input_data_dir + file_name).c_str(), "r");

            while (fscanf(level_scs[i], "%d %d %d", &u, &v, &dis) == 3) {

                int flg = 0;
                for (KNode knode : level_scs_arr[i][u]) {
                    if (knode.id == v && knode.dis == dis) {
                        flg = 1;
                        break;
                    }
                }
                if (flg == 0) {
                    if (u > test_n) {
                        cout << "during reading: u>test_n u: " << u << endl;
                        stop_here();
                    }
                    if (v > test_n) {
                        cout << "during reading: v>test_n v: " << v << endl;
                        stop_here();
                    }
                    KNode *newKnode = new KNode(v, dis);
                    level_scs_arr[i][u].push_back(*newKnode);
                }

            }
            fclose(level_scs[i]);
        }

        for (int i = 0; i <= lv_num; i++) {
            for (int j = 0; j < MAGIC_NUM; j++) {
                vector<KNode> nei = level_scs_arr[i][j];
                for (KNode knode : nei) {
                    int v = knode.id;
                    int d = knode.dis;
                    if (v > test_n) {
                        cout << "during testing: " << i << " " << j << " " << node_order[j] << " " << v << " " << d
                             << endl;
                        stop_here();
                    }
                }
            }
        }
        for (int i = 0; i <= test_n; i++) {
            for (KNode knode : rev_shortcuts_arr[i]) {
                int v = knode.id;
                int d = knode.dis;
                if (v > test_n) {
                    cout << "read rev_shortcuts_arr error: " << i << endl;
                    stop_here();
                }
            }
        }
        for (int i = 0; i <= test_n; i++) {
            for (KNode knode : rev_shortcuts_level_arr[i]) {
                int v = knode.id;
                int d = knode.dis;
                if (v > test_n) {
                    cout << "read rev_shortcuts_level_arr error: " << i << endl;
                    stop_here();
                }
            }
        }
        cout << "read shortcut ends..." << endl;
        has_read_scob_index=1;
    }
}



void update_hier_local_knn_carpos_edgeversion(CarPos &cp, int k_max, long *dis, int *visited, DijkstraQueue *q) {
    vector<int> visitedVector;

    dis[cp.node1] = cp.dis1;
    q->push(cp.node1, cp.dis1);
    visited[cp.node1] = 1;
    visitedVector.push_back(cp.node1);

    dis[cp.node2] = cp.dis2;
    q->push(cp.node2, cp.dis2);
    visited[cp.node2] = 1;
    visitedVector.push_back(cp.node2);


    while (!q->empty()) {
        int u = q->pop();
        if (u > test_n)
            continue;
        vector<KNode> knodes = hier_local_knn_arr[u];
        int j;
        for (j = 0; j < knodes.size(); j++) {
            if (dis[u] < knodes[j].dis)
                break;
        }
//		for (j = knodes.size()-1; j >=0; j--) {
//					if (dis[u] >= knodes[j].dis)
//						break;
//				}
//		j++;
        if (j < k_max) {
            if (knodes.size() < k_max) {
                KNode *c = new KNode(test_n + cp.id, dis[u]);
                knodes.push_back(*c);
            }

            int size = knodes.size();
            for (int t = size - 1; t > j; t--)
                knodes[t] = knodes[t - 1];
            knodes[j].id = test_n + cp.id;
            knodes[j].dis = dis[u];
//			car_localarr_map[cp.id].push_back(u);
            hier_local_knn_arr[u] = knodes;
        } else {
            hier_local_knn_arr[u] = knodes;
            continue;

        }

//		hier_local_knn_arr[u] = knodes;
        int need_expand = 1;
//		cout<<"here1"<<endl;
        vector<KNode> nei;
        if (mode == 1) {
            if (node_order[u] >= cut_level) {
                continue;
            } else if (node_order[u] >= tradeoff_level)
                nei = shortcuts_inc_arr[u];
            else
                nei = shortcuts_arr[u];
        } else if (mode == 2) {
            if (node_order[u] >= cut_level) {
                nei = level_scs_arr[cut_level][u];
                //				cout<<nei.size()<<endl;
            } else if (node_order[u] >= tradeoff_level)
                nei = shortcuts_arr[u];
            else
                nei = shortcuts_inc_arr[u];
        } else {
            if (node_order[u] >= cut_level) {
                continue;
            } else if (node_order[u] >= tradeoff_level)
                nei = shortcuts_arr[u];
            else
                nei = shortcuts_inc_arr[u];
        }
        if (cut_level) {
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
//				int level_gap_for_v = level_gap[node_order[v]] * 2 + 1;
//				int coor_v_x = cell_loc_x[v];
//				int coor_v_y = cell_loc_y[v];
//				int gap1 = max(abs(coor_x1 - coor_v_x),
//						abs(coor_y1 - coor_v_y));
//				int gap2 = max(abs(coor_x2 - coor_v_x),
//						abs(coor_y2 - coor_v_y));
//				if (gap1 > level_gap_for_v && gap2 > level_gap_for_v)
//					continue;
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
    for (int item : visitedVector) {
        visited[item] = 0;
    }

    q->clear();
}

void update_hier_local_knn_carpos_nodeversion_enu(int is_init, int nodeid, int k_max, int k) {


    for (KNode knode: shortcuts_inc_arr_local[cut_level][nodeid]) {
        int u = knode.id;
        int dis_u = knode.dis;
        vector<KNode> &knodes = hier_local_knn_arr[u];
        int j;
        for (j = 0; j < knodes.size(); j++) {
            if (dis_u < knodes[j].dis)
                break;
        }
        if (is_init) {
            if (j < k_max) {
                if (knodes.size() < k_max) {
                    KNode *c = new KNode(nodeid, dis_u);
                    knodes.push_back(*c);
                }

                int size = knodes.size();
                for (int t = size - 1; t > j; t--)
                    knodes[t] = knodes[t - 1];
                knodes[j].id = nodeid;
                knodes[j].dis = dis_u;
            } else {

                continue;

            }
        }
        if (!is_init) {
            if (j < knodes.size() || (j < k_max && knodes.size() < k)) {
                // if knodes.size()<k, then the list is correct
                if (knodes.size() < k_max) {
                    KNode *c = new KNode(nodeid, dis_u);
                    knodes.push_back(*c);
                }
                int size = knodes.size();
                for (int t = size - 1; t > j; t--)
                    knodes[t] = knodes[t - 1];
                knodes[j].id = nodeid;
                knodes[j].dis = dis_u;
            }

        }


    }

}


void update_kDNN_object(int is_init, int nodeid, int k_max, int k, long *dis, int *visited,
                        DijkstraQueue *q) {
//	cout<<"start update..."<<endl;
    vector<int> visitedVector;
//	int coor_x = cell_loc_x[nodeid];
//	int coor_y = cell_loc_y[nodeid];

    dis[nodeid] = 0;
    q->push(nodeid, 0);
    visited[nodeid] = 1;
    visitedVector.push_back(nodeid);

    while (!q->empty()) {
        int u = q->pop();
        vector<KNode> &knodes = hier_local_knn_arr[u];
        int j;
        for (j = 0; j < knodes.size(); j++) {
            if (dis[u] < knodes[j].dis)
                break;
        }
        if (is_init) {
            if (j < k_max) {
                if (knodes.size() < k_max) {
                    KNode *c = new KNode(nodeid, dis[u]);
                    knodes.push_back(*c);
                }

                int size = knodes.size();
                for (int t = size - 1; t > j; t--)
                    knodes[t] = knodes[t - 1];
                knodes[j].id = nodeid;
                knodes[j].dis = dis[u];
            } else {

                continue;

            }
        }
        if (!is_init) {
            // further optimization is possible here
            if (j < knodes.size() || (j < k_max && knodes.size() < k)) {
                // when insert, not to exceed a length of k_max
                // if knodes.size()<k, then the list is correct
                if (knodes.size() < k_max) {
                    KNode *c = new KNode(nodeid, dis[u]);
                    knodes.push_back(*c);
                }
                int size = knodes.size();
                for (int t = size - 1; t > j; t--)
                    knodes[t] = knodes[t - 1];
                knodes[j].id = nodeid;
                knodes[j].dis = dis[u];
            }

        }
        int need_expand = 1;

        vector<KNode> *nei;

        if (mode == 2) {
            if (node_order[u] >= cut_level) {
                nei = &(level_scs_arr[cut_level][u]);
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
            else
                nei = &(shortcuts_inc_arr[u]);
//				nei = &(shortcuts_inc_arr_local[cut_level][u]);

        } else {
            if (node_order[u] >= cut_level) {
                continue;
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
            else
                nei = &(shortcuts_inc_arr[u]);
//				nei = &(shortcuts_inc_arr_local[cut_level][u]);

        }

        if (cut_level) {
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

//				int level_gap_for_v = level_gap[node_order[v]] * 2 + 1;
//				int coor_v_x = cell_loc_x[v];
//				int coor_v_y = cell_loc_y[v];
//				int gap = max(abs(coor_x - coor_v_x), abs(coor_y - coor_v_y));

//				if (gap > level_gap_for_v)
//					continue;
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
    for (int item : visitedVector) {
        visited[item] = 0;
    }

    q->clear();
}

void update_kDNN_object_multi(int is_init, int nodeid, int k_max, int k, long *dis, int *visited,
                        DijkstraQueue *q, vector<KNode>* hier_local_knn_arr) {
//	cout<<"start update..."<<endl;
    vector<int> visitedVector;
//	int coor_x = cell_loc_x[nodeid];
//	int coor_y = cell_loc_y[nodeid];

    dis[nodeid] = 0;
    q->push(nodeid, 0);
    visited[nodeid] = 1;
    visitedVector.push_back(nodeid);
//    cout<<"start exploration..."<<endl;
    while (!q->empty()) {
//        cout<<"here 1"<<endl;
        int u = q->pop();
//        cout<<"u: "<<u<<endl;
        vector<KNode> &knodes = hier_local_knn_arr[u];
        int j=0;
        for (j = 0; j < knodes.size(); j++) {
            if (dis[u] < knodes[j].dis)
                break;
        }
//        cout<<"here 2"<<endl;
        if (is_init) {
            if (j < k_max) {
                if (knodes.size() < k_max) {
//                    cout<<" before new "<<endl;
                    KNode *c = new KNode(nodeid, dis[u]);
//                    cout<<"after new"<<endl;
                    knodes.push_back(*c);
//                    cout<<"after push back"<<endl;
                }

                int size = knodes.size();
                for (int t = size - 1; t > j; t--)
                    knodes[t] = knodes[t - 1];
                knodes[j].id = nodeid;
                knodes[j].dis = dis[u];
            } else {

                continue;

            }
        }
        if (!is_init) {
            // further optimization is possible here
//            cout<<"k_max: "<<k_max<<endl;
//            cout<<"knodes.size(): ";
//            cout<<knodes.size()<<endl;
//            cout<<"j: "<<endl;
//            cout<<j<<endl;
            if (j < knodes.size() || (j < k_max && knodes.size() < k)) {
                // when insert, not to exceed a length of k_max
                // if knodes.size()<k, then the list is correct
                if (knodes.size() < k_max) {
//                    cout<<"nodeid, dis[u]: "<<nodeid<<" "<<dis[u]<<endl;
                    KNode *c = new KNode(nodeid, dis[u]);
//                    cout<<"after new"<<endl;
                    knodes.push_back(*c);
//                    cout<<"after push back"<<endl;
                }
                int size = knodes.size();
//                cout<<"size: "<<size<<endl;
                for (int t = size - 1; t > j; t--)
                    knodes[t] = knodes[t - 1];
                knodes[j].id = nodeid;
                knodes[j].dis = dis[u];
            }

        }
//        cout<<"here 3"<<endl;
        int need_expand = 1;

        vector<KNode> *nei;

        if (mode == 2) {
//            cout<<"there 1"<<endl;
//            cout<<"u: "<<endl;
//            cout<<u<<endl;
//            cout<<"node_order[u]: "<<endl;
//            cout<<node_order[u]<<endl;


            if (node_order[u] >= cut_level) {
//                cout<<"level_scs_arr[cut_level][u] size: "<<endl;
//                cout<<level_scs_arr[cut_level][u].size()<<endl;
                nei = &(level_scs_arr[cut_level][u]);
            } else if (node_order[u] >= tradeoff_level) {
//                cout<<"shortcuts_arr[u]: "<<endl;
//                cout<<shortcuts_arr[u].size()<<endl;
                nei = &(shortcuts_arr[u]);
            }
            else {
//                cout<<"shortcuts_inc_arr[u]: "<<endl;
//                cout<<shortcuts_inc_arr[u].size()<<endl;
                nei = &(shortcuts_inc_arr[u]);
//				nei = &(shortcuts_inc_arr_local[cut_level][u]);
            }

        } else {
            if (node_order[u] >= cut_level) {
                continue;
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
            else
                nei = &(shortcuts_inc_arr[u]);
//				nei = &(shortcuts_inc_arr_local[cut_level][u]);

        }
//        cout<<"here 4"<<endl;
//        cout<<"size: "<<endl;
//        cout<<nei->size()<<endl;

        if (cut_level) {
            for (KNode knode : *nei) {
                int v = knode.id;
                if (visited[v] && dis[v] + knode.dis < dis[u]) {
                    dis[u] = dis[v] + knode.dis;
                    need_expand = 0;
                    break;
                }
            }
        }

//        cout<<"here 5"<<endl;

        if (need_expand) {
            for (KNode knode : *nei) {

                int v = knode.id;

//				int level_gap_for_v = level_gap[node_order[v]] * 2 + 1;
//				int coor_v_x = cell_loc_x[v];
//				int coor_v_y = cell_loc_y[v];
//				int gap = max(abs(coor_x - coor_v_x), abs(coor_y - coor_v_y));

//				if (gap > level_gap_for_v)
//					continue;
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
//        cout<<"here 6"<<endl;

    }
//    cout<<"end exploration..."<<endl;
    for (int item : visitedVector) {
        visited[item] = 0;
    }

    q->clear();
}

vector<long> update_kdnn_list(int is_init, int k_max, int k, int version, DijkstraQueue *q, long *dis,
                              int *visited, vector<int> &node_ran_objs) {
    vector<long> cgraph_each_update_time;

    memset(visited, 0, sizeof(int) * (test_n + 1 + test_n));
    for (int nodeid : node_ran_objs) {
        long start = clock();
        update_kDNN_object(is_init, nodeid, k_max, k, dis, visited, q);
        cgraph_each_update_time.push_back(clock() - start);
    }
    return cgraph_each_update_time;

}

void update_dijkstra(unordered_map<int, vector<PointCar> > &car_hash) {
    car_hash.clear();
    for (CarPos cp : carposes) {
        if (car_hash.find(cp.node1) == car_hash.end()) {
            vector<PointCar> temp;
            car_hash[cp.node1] = temp;
        }
        car_hash[cp.node1].push_back(PointCar(cp.id, cp.dis1));

        if (car_hash.find(cp.node2) == car_hash.end()) {
            vector<PointCar> temp;
            car_hash[cp.node2] = temp;
        }
        car_hash[cp.node2].push_back(PointCar(cp.id, cp.dis2));

    }

}


vector<KNode> updateToList(vector<KNode> &knodes, int carid, int dis, int k_max) {
    int i = 0;
    for (i = 0; i < knodes.size(); i++) {
        if (knodes[i].id == carid)
            break;
    }
    if (i == knodes.size()) {
        //not exists

        int j;
        for (j = 0; j < knodes.size(); j++) {
            if (dis < knodes[j].dis) {
                break;
            }
        }
        if (j < k_max) {
            if (knodes.size() < k_max) {
                KNode *c = new KNode(carid, dis);
                knodes.push_back(*c);
            }

            int size = knodes.size();
            for (int t = size - 1; t > j; t--) {
                knodes[t] = knodes[t - 1];
            }
            knodes[j].id = carid;
            knodes[j].dis = dis;

        }

    } else {
        if (knodes[i].dis > dis) {
            knodes[i].dis = dis;
            for (int j = i; j > 0; j--) {
                if (knodes[j].dis < knodes[j - 1].dis) {
                    KNode tmp = knodes[j];
                    knodes[j] = knodes[j - 1];
                    knodes[j - 1] = tmp;
                }
            }
        }
    }
    return knodes;
}

vector<KNode> mergers(vector<KNode> &knodes1, vector<KNode> &knodes2, int k_max, int dis) {
    vector<KNode> rs;
    int i = 0;
    int j = 0;
    int distinct = 0;
    unordered_map<int, int> hash;
    while (distinct < k_max && (i < knodes1.size() || j < knodes2.size())) {
        if (i == knodes1.size() && j < knodes2.size()) {

            if (hash.find(knodes2[j].id) == hash.end()) {
                hash[knodes2[j].id] = 1;
                rs.push_back(KNode(knodes2[j].id, knodes2[j].dis + dis));
                distinct++;
            }
            j++;

        }
        if (i < knodes1.size() && j == knodes2.size()) {

            if (hash.find(knodes1[i].id) == hash.end()) {
                hash[knodes1[i].id] = 1;
                rs.push_back(knodes1[i]);

                distinct++;
            }
            i++;
        }
        if (i < knodes1.size() && j < knodes2.size()) {
            if (knodes1[i].dis <= knodes2[j].dis + dis) {
                if (hash.find(knodes1[i].id) == hash.end()) {
                    rs.push_back(knodes1[i]);
                    hash[knodes1[i].id] = 1;

                    distinct++;
                }
                i++;
            }
            if (knodes1[i].dis > knodes2[j].dis + dis) {
                if (hash.find(knodes2[j].id) == hash.end()) {
                    rs.push_back(KNode(knodes2[j].id, knodes2[j].dis + dis));
                    hash[knodes2[j].id] = 1;

                    distinct++;
                }
                j++;
            }
        }
    }
    return rs;
}

vector<KNode> update_node_downward(int k, int query, int avoidid) {

    vector<KNode> down_nei = rev_shortcuts_arr[query];
    vector<KNode> list = hier_local_knn_arr[query];
    for (KNode knode : down_nei) {
        if (node_order[query] < node_order[knode.id]) {
            cout << "error1" << endl;
        }
        if (node_order[query] == node_order[knode.id] && inlev_order[query] < inlev_order[knode.id]) {
            cout << "error2" << endl;
        }
        if (node_order[query] == node_order[knode.id] && inlev_order[query] == inlev_order[knode.id]
            && query < knode.id) {
            cout << "error3" << endl;
        }
        for (KNode res : hier_local_knn_arr[knode.id]) {
//			if (res.id == avoidid) {
//				cout << "bb " << avoidid - test_n << endl;
//				cout << "node: " << knode.id << endl;
//				stop_here();
//			}
            list = updateToList(list, res.id, res.dis + knode.dis, k);
        }
    }
    return list;
}

void update_node_downward_node(int k, int query, DijkstraQueue *q, DijkstraQueue *q1, long *dis, long *dis1,
                               int *visited, int *colors, int *marked_first, vector<KNode> *&rev_shortcuts_arr) {

    //initially, if(marked_first[query]==0), then the node $query$'s list is correct

    /*
     * Pruning 1: We use variable $colors for an efficient stopping rule.
     * (1) Initially all the nodes are colored "0" except the source node $query$.
     * (2) If the (downward) search started at $query$ now visits a node $v$ that has a different list from that of $query$, then
     * we know that the top-k candidates through the node $v$ must be found.
     * (3) If the search boundaries are all colored "0", we can safely stop the downward search started at u.
     */
    int rank = node_order[query];

    // preparing for downward search started at $query$
    vector<int> visitedVector;
    dis[query] = 0;
    q->push(query, 0);
    visited[query] = 1;
    visitedVector.push_back(query);

    // q1 maintains a distance lower bound of the kNN of $query$
    // note: the size of $result is k-1
    vector<KNode> &result = hier_local_knn_arr[query];
    if (result.size()) {
        q1->push(query, result[result.size() - 1].dis);
        dis1[query] = result[result.size() - 1].dis;
    } else {
        q1->push(query, dis[query]);
        dis1[query] = dis[query];
    }

    // use $sum$ as the code to compare two lists
    int sum = 0;
    int query_list_size = result.size();
    for (KNode knode : result)
        sum += knode.id;
    int unmarked_cnt = 1;
    int flag_boundary_meet = 0;
    colors[query] = 1;
    // dis1[min_k_dis_boundary] is a reasonable lower bound of the searched radius
    int min_k_dis_boundary = query;
    while (!q->empty()) {
        int u = q->pop();
        // note that, we use a DijkstraQueue q1 (not q!!) to maintain a lowerbound of the top-k object distance of the query node
        if (flag_boundary_meet == 0) {
            while (!q1->empty()) {
                min_k_dis_boundary = q1->top();
                if (visited[min_k_dis_boundary] == 2) {
                    q1->pop();
                } else {
                    break;
                }
            }
            visited[u] = 2;
        }
        int need_expand = 1;
        vector<KNode> &knodes = hier_local_knn_arr[u];
        int size = knodes.size();

        // use the list of u to update the final result
        // only maintain the correctness of top-k
        int t = 0;
        for (int i = 0; i < k && i < size; ++i) {
            KNode &knode = knodes[i];
            int updated_dis = dis[u] + knode.dis;
            if (result.size() >= k && updated_dis >= result[k - 1].dis) {
                break;
            }
            for (; t < result.size(); t++) {
                if (updated_dis < result[t].dis)
                    break;
            }

            // $t$ now stores the index of $result$ where the new value $updated_dis$ should be inserted.

            if (t == k) {
                break;
            }
            int t2 = 0;
            for (t2 = result.size() - 1; t2 >= 0; t2--) {
                if (result[t2].id == knode.id)
                    break;
            }

            // $t2$ now stores the index in $result$ where the object id is the same as the one to be inserted. That means
            // we need to update the distance value of that object using $updated_dis$
            if (t2 == -1)
                t2 = result.size();
            if (t2 == result.size()) {
                // if the object is new to $result$
                if (result.size() < k) {
                    KNode *c = new KNode(knode.id, updated_dis);
                    result.push_back(*c);
                }
                for (int h = result.size() - 1; h > t; h--)
                    result[h] = result[h - 1];
                if (t < k) {
                    result[t].id = knode.id;
                    result[t].dis = updated_dis;
                }

            } else {
                // if the object has already been in $result$
//                int t1 = max(t, t2);
//                for (; t1 < result.size(); t1++) {
//                    if (result[t1].id == knode.id && updated_dis < result[t1].dis)
//                        break;
//                }
//                if (t1 < result.size()) {
//                    for (int h = t1; h > t; h--)
//                        result[h] = result[h - 1];
//                    if (t < k) {
//                        result[t].id = knode.id;
//                        result[t].dis = updated_dis;
//
//                    }
//                }

                if(t<=t2){
                    for (int h = t2; h > t; h--)
                        result[h] = result[h - 1];
                    if (t < k) {
                        result[t].id = knode.id;
                        result[t].dis = updated_dis;
                    }

                }
            }

        }
        // Pruning rule 1: unmarked_cnt denotes the number of nodes in the search boundaries that are colored "1"
        if (unmarked_cnt == 0)
            // if the search is enclosed by all nodes with colors of 0, then go on the search is meaningless, but
//            still needs to explore the maintained boundaries.
            continue;
        if (colors[u] == 1) {
            // pop out node u, and reduce the unmarked_cnt if colors[u]=1
            unmarked_cnt--;
        }
        // Pruning rule 2: if
        if (flag_boundary_meet)
            continue;
        if (result.size() >= k && dis1[min_k_dis_boundary] > result[k - 1].dis) {
            // If goes on expansion, dis1[min_k_dis_boundary] is the best distance
            flag_boundary_meet = 1;
            continue;
        }

        vector<KNode> &nei = rev_shortcuts_arr[u];
//        vector<KNode> & nei = rev_shortcuts_level_arr[u];

//		vector<KNode>* nei=NULL;
//		vector<KNode> list;
//
//		if (mode == 2) {
//			if (rank >= cut_level) {
//				for (KNode knode : level_scs_arr[cut_level][u])
//					list.push_back(knode);
//				for (KNode knode : rev_shortcuts_level_arr[u])
//					if (node_order[knode.id] < cut_level)
//						list.push_back(knode);
//				for (KNode knode : rev_shortcuts_arr[u]) {
//					if (node_order[knode.id] < tradeoff_level)
//						list.push_back(knode);
//				}
//				nei = &list;
//			} else if (rank >= tradeoff_level) {
//				for (KNode knode : rev_shortcuts_level_arr[u])
//					list.push_back(knode);			// level decrease
//				for (KNode knode : rev_shortcuts_arr[u]) {
//					if (node_order[knode.id] < tradeoff_level)
//						list.push_back(knode);			// strict decrease
//				}
//				nei = &list;
//			} else
//				nei = &(rev_shortcuts_arr[u]);			//strict decrease
//		} else {
//			if (rank >= cut_level) {
//				for (KNode knode : rev_shortcuts_level_arr[u]) {
//					if (node_order[knode.id] < cut_level)
//						list.push_back(knode);
//				}
//				for (KNode knode : rev_shortcuts_arr[u])
//					if (node_order[knode.id] < tradeoff_level)
//						list.push_back(knode);
//				nei = &list;
//
//			} else if (rank >= tradeoff_level) {
//				for (KNode knode : rev_shortcuts_level_arr[u])
//					list.push_back(knode);
//				for (KNode knode : rev_shortcuts_arr[u])
//					if (node_order[knode.id] < tradeoff_level)
//						list.push_back(knode);
//				nei = &list;
//			} else
//				nei = &(rev_shortcuts_arr[u]);
//		}

        if (cut_level) {
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
                vector<KNode> &li = hier_local_knn_arr[v];

                if (visited[v]) {
                    if (knode.dis < dis[v] - dis[u]) {
                        dis[v] = dis[u] + knode.dis;

                        // to maintain top-k object distance of query node (note: using q1, a different DijkstraQueue)
                        int entry = li.size() - 1;
                        if (entry > k - 2)
                            entry = k - 1;
                        if (entry >= 0) {
                            dis1[v] = dis[v] + li[entry].dis;
                            q1->decrease_priority(v, dis1[v]);

                        } else {
                            dis1[v] = dis[v];
                            q1->decrease_priority(v, dis1[v]);
                        }
                        q->decrease_priority(v, dis[v]);
                        //color's algorithm to test whether the boundaries of the search range are all stable nodes, that is, with correct object
                        // lists.
                        if (colors[v])
                            unmarked_cnt--;
                        colors[v] = colors[u];
                        if (colors[v])
                            unmarked_cnt++;
                    }
                } else {
                    dis[v] = dis[u] + knode.dis;
                    q->push(v, dis[v]);
                    int entry = li.size() - 1;
                    if (entry > k - 2)
                        entry = k - 1;
                    if (entry >= 0) {
                        dis1[v] = dis[v] + li[entry].dis;
                        q1->push(v, dis1[v]);

                    } else {
                        dis1[v] = dis[v];
                        q1->push(v, dis1[v]);
                    }

                    visited[v] = 1;
                    visitedVector.push_back(v);
                    int flag = 1;
                    if (marked_first[v] == 0) {
                        // It means the node v has a correct associated list
                        flag = 0;
                    }
                    else if (k - 1 > li.size()) {
                        flag = 0; // It means the node v has a correct top-k list
                    }
                    else if (query_list_size == li.size()) {
                        int sum_v = 0;
                        for (KNode knode : hier_local_knn_arr[v])
                            sum_v += knode.id;
                        if (sum != sum_v)
                            // It means the list of v is different from the list of query, since the two sizes of two lists
                            // are the same, while the lists have different ID sums
                            flag = 0;
                    }
                    if (flag) {
                        colors[v] = 1;
                        unmarked_cnt++;
                    }

                }

            }

        }

    }
    for (int item : visitedVector) {
        visited[item] = 0;
        colors[item] = 0;
    }
    q->clear();
    q1->clear();
}

void update_node_downward_node_multi(int k, int query, DijkstraQueue *q, DijkstraQueue *q1, long *dis, long *dis1,
                               int *visited, int *colors, int *marked_first, vector<KNode> *&rev_shortcuts_arr, vector<KNode>* hier_local_knn_arr) {

    //initially, if(marked_first[query]==0), then the node $query$'s list is correct

    /*
     * Pruning 1: We use variable $colors for an efficient stopping rule.
     * (1) Initially all the nodes are colored "0" except the source node $query$.
     * (2) If the (downward) search started at $query$ now visits a node $v$ that has a different list from that of $query$, then
     * we know that the top-k candidates through the node $v$ must be found.
     * (3) If the search boundaries are all colored "0", we can safely stop the downward search started at u.
     */
    int rank = node_order[query];

    // preparing for downward search started at $query$
    vector<int> visitedVector;
    dis[query] = 0;
    q->push(query, 0);
    visited[query] = 1;
    visitedVector.push_back(query);

    // q1 maintains a distance lower bound of the kNN of $query$
    // note: the size of $result is k-1
    vector<KNode> &result = hier_local_knn_arr[query];
    if (result.size()) {
        q1->push(query, result[result.size() - 1].dis);
        dis1[query] = result[result.size() - 1].dis;
    } else {
        q1->push(query, dis[query]);
        dis1[query] = dis[query];
    }

    // use $sum$ as the code to compare two lists
    int sum = 0;
    int query_list_size = result.size();
    for (KNode knode : result)
        sum += knode.id;
    int unmarked_cnt = 1;
    int flag_boundary_meet = 0;
    colors[query] = 1;
    // dis1[min_k_dis_boundary] is a reasonable lower bound of the searched radius
    int min_k_dis_boundary = query;
    while (!q->empty()) {
        int u = q->pop();
        // note that, we use a DijkstraQueue q1 (not q!!) to maintain a lowerbound of the top-k object distance of the query node
        if (flag_boundary_meet == 0) {
            while (!q1->empty()) {
                min_k_dis_boundary = q1->top();
                if (visited[min_k_dis_boundary] == 2) {
                    q1->pop();
                } else {
                    break;
                }
            }
            visited[u] = 2;
        }
        int need_expand = 1;
        vector<KNode> &knodes = hier_local_knn_arr[u];
        int size = knodes.size();

        // use the list of u to update the final result
        // only maintain the correctness of top-k
        int t = 0;
        for (int i = 0; i < k && i < size; ++i) {
            KNode &knode = knodes[i];
            int updated_dis = dis[u] + knode.dis;
            if (result.size() >= k && updated_dis >= result[k - 1].dis) {
                break;
            }
            for (; t < result.size(); t++) {
                if (updated_dis < result[t].dis)
                    break;
            }

            // $t$ now stores the index of $result$ where the new value $updated_dis$ should be inserted.

            if (t == k) {
                break;
            }
            int t2 = 0;
            for (t2 = result.size() - 1; t2 >= 0; t2--) {
                if (result[t2].id == knode.id)
                    break;
            }

            // $t2$ now stores the index in $result$ where the object id is the same as the one to be inserted. That means
            // we need to update the distance value of that object using $updated_dis$
            if (t2 == -1)
                t2 = result.size();
            if (t2 == result.size()) {
                // if the object is new to $result$
                if (result.size() < k) {
                    KNode *c = new KNode(knode.id, updated_dis);
                    result.push_back(*c);
                }
                for (int h = result.size() - 1; h > t; h--)
                    result[h] = result[h - 1];
                if (t < k) {
                    result[t].id = knode.id;
                    result[t].dis = updated_dis;
                }

            } else {
                // if the object has already been in $result$
//                int t1 = max(t, t2);
//                for (; t1 < result.size(); t1++) {
//                    if (result[t1].id == knode.id && updated_dis < result[t1].dis)
//                        break;
//                }
//                if (t1 < result.size()) {
//                    for (int h = t1; h > t; h--)
//                        result[h] = result[h - 1];
//                    if (t < k) {
//                        result[t].id = knode.id;
//                        result[t].dis = updated_dis;
//
//                    }
//                }

                if(t<=t2){
                    for (int h = t2; h > t; h--)
                        result[h] = result[h - 1];
                    if (t < k) {
                        result[t].id = knode.id;
                        result[t].dis = updated_dis;
                    }

                }
            }

        }
        // Pruning rule 1: unmarked_cnt denotes the number of nodes in the search boundaries that are colored "1"
        if (unmarked_cnt == 0)
            // if the search is enclosed by all nodes with colors of 0, then go on the search is meaningless, but
//            still needs to explore the maintained boundaries.
            continue;
        if (colors[u] == 1) {
            // pop out node u, and reduce the unmarked_cnt if colors[u]=1
            unmarked_cnt--;
        }
        // Pruning rule 2: if
        if (flag_boundary_meet)
            continue;
        if (result.size() >= k && dis1[min_k_dis_boundary] > result[k - 1].dis) {
            // If goes on expansion, dis1[min_k_dis_boundary] is the best distance
            flag_boundary_meet = 1;
            continue;
        }

        vector<KNode> &nei = rev_shortcuts_arr[u];
//        vector<KNode> & nei = rev_shortcuts_level_arr[u];

//		vector<KNode>* nei=NULL;
//		vector<KNode> list;
//
//		if (mode == 2) {
//			if (rank >= cut_level) {
//				for (KNode knode : level_scs_arr[cut_level][u])
//					list.push_back(knode);
//				for (KNode knode : rev_shortcuts_level_arr[u])
//					if (node_order[knode.id] < cut_level)
//						list.push_back(knode);
//				for (KNode knode : rev_shortcuts_arr[u]) {
//					if (node_order[knode.id] < tradeoff_level)
//						list.push_back(knode);
//				}
//				nei = &list;
//			} else if (rank >= tradeoff_level) {
//				for (KNode knode : rev_shortcuts_level_arr[u])
//					list.push_back(knode);			// level decrease
//				for (KNode knode : rev_shortcuts_arr[u]) {
//					if (node_order[knode.id] < tradeoff_level)
//						list.push_back(knode);			// strict decrease
//				}
//				nei = &list;
//			} else
//				nei = &(rev_shortcuts_arr[u]);			//strict decrease
//		} else {
//			if (rank >= cut_level) {
//				for (KNode knode : rev_shortcuts_level_arr[u]) {
//					if (node_order[knode.id] < cut_level)
//						list.push_back(knode);
//				}
//				for (KNode knode : rev_shortcuts_arr[u])
//					if (node_order[knode.id] < tradeoff_level)
//						list.push_back(knode);
//				nei = &list;
//
//			} else if (rank >= tradeoff_level) {
//				for (KNode knode : rev_shortcuts_level_arr[u])
//					list.push_back(knode);
//				for (KNode knode : rev_shortcuts_arr[u])
//					if (node_order[knode.id] < tradeoff_level)
//						list.push_back(knode);
//				nei = &list;
//			} else
//				nei = &(rev_shortcuts_arr[u]);
//		}

        if (cut_level) {
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
                vector<KNode> &li = hier_local_knn_arr[v];

                if (visited[v]) {
                    if (knode.dis < dis[v] - dis[u]) {
                        dis[v] = dis[u] + knode.dis;

                        // to maintain top-k object distance of query node (note: using q1, a different DijkstraQueue)
                        int entry = li.size() - 1;
                        if (entry > k - 2)
                            entry = k - 1;
                        if (entry >= 0) {
                            dis1[v] = dis[v] + li[entry].dis;
                            q1->decrease_priority(v, dis1[v]);

                        } else {
                            dis1[v] = dis[v];
                            q1->decrease_priority(v, dis1[v]);
                        }
                        q->decrease_priority(v, dis[v]);
                        //color's algorithm to test whether the boundaries of the search range are all stable nodes, that is, with correct object
                        // lists.
                        if (colors[v])
                            unmarked_cnt--;
                        colors[v] = colors[u];
                        if (colors[v])
                            unmarked_cnt++;
                    }
                } else {
                    dis[v] = dis[u] + knode.dis;
                    q->push(v, dis[v]);
                    int entry = li.size() - 1;
                    if (entry > k - 2)
                        entry = k - 1;
                    if (entry >= 0) {
                        dis1[v] = dis[v] + li[entry].dis;
                        q1->push(v, dis1[v]);

                    } else {
                        dis1[v] = dis[v];
                        q1->push(v, dis1[v]);
                    }

                    visited[v] = 1;
                    visitedVector.push_back(v);
                    int flag = 1;
                    if (marked_first[v] == 0) {
                        // It means the node v has a correct associated list
                        flag = 0;
                    }
                    else if (k - 1 > li.size()) {
                        flag = 0; // It means the node v has a correct top-k list
                    }
                    else if (query_list_size == li.size()) {
                        int sum_v = 0;
                        for (KNode knode : hier_local_knn_arr[v])
                            sum_v += knode.id;
                        if (sum != sum_v)
                            // It means the list of v is different from the list of query, since the two sizes of two lists
                            // are the same, while the lists have different ID sums
                            flag = 0;
                    }
                    if (flag) {
                        colors[v] = 1;
                        unmarked_cnt++;
                    }

                }

            }

        }

    }
    for (int item : visitedVector) {
        visited[item] = 0;
        colors[item] = 0;
    }
    q->clear();
    q1->clear();
}


int remove_kdnn_object(int remove_node, int k_max, int k, mem_struct &mems,
                       vector<KNode> *&rev_shortcuts_arr,
                       int durationLimits) {
    int durationLimitsMicro = durationLimits * 1000000;
    vector<int> visitedVector;
    DijkstraQueue *q = mems.q;
    DijkstraQueue *q1 = mems.q1;
    long *dis = mems.dist;
    long *dis1 = mems.dist1;
    int *visited = mems.visited;
    int *marked_first = mems.marked_first;
    int *colors = mems.visited1;
    dis[remove_node] = 0;
    q->push(remove_node, 0);
    visited[remove_node] = 1;
    visitedVector.push_back(remove_node);

    vector<NodeRank> nodecoll;
    while (!q->empty()) {
        int u = q->pop();
        if (u > test_n)
            continue;

        vector<KNode> &knodes = hier_local_knn_arr[u];
        vector<KNode>::iterator it;
        int is_erase = 0;
        for (it = knodes.begin(); it != knodes.end(); it++) {
            if (remove_node == it->id) {
                it = knodes.erase(it);
                is_erase = 1;
                break;
            }
        }
        // If originally, the size of the object list is less than k, that means, after object removal, it remains
        // less than k-1 objects, then all the downhill objects are already in the list. Hence, we only push in the node
        // when knodes.size()==k-1
        // set marked_first[u]=1, means that the current list is unstable, due to the object removal
        if (is_erase && knodes.size() == k - 1) {
            nodecoll.push_back(NodeRank(node_order[u], inlev_order[u], u));
            marked_first[u] = 1;
        }


        vector<KNode> *nei;
        if (mode == QUERY_FAVOR_MODE) {
            if (node_order[u] >= cut_level) {
                nei = &(level_scs_arr[cut_level][u]);
                /*
                 * if node u'level >= cut_level, then associate the node u with a level $cut_level, and employ its level
                 * shortcuts.
                 */
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
                /*
                 * if at least tradeoff_level, then "level+up"
                 */
            else
                nei = &(shortcuts_inc_arr[u]);
            /*
             * if less that tradeoff_level, then "up"
             */
        } else {
            if (node_order[u] >= cut_level) {
                continue;
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
            else
                nei = &(shortcuts_inc_arr[u]);
        }


        for (KNode knode : *nei) {
            int v = knode.id;
//				int level_gap_for_v = level_gap[node_order[v]] * 2 + 1;
//				int coor_v_x = cell_loc_x[v];
//				int coor_v_y = cell_loc_y[v];
//				int gap1 = max(abs(coor_x1 - coor_v_x),
//						abs(coor_y1 - coor_v_y));
//				int gap2 = max(abs(coor_x2 - coor_v_x),
//						abs(coor_y2 - coor_v_y));
//				if (gap1 > level_gap_for_v && gap2 > level_gap_for_v)
//					continue;
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
    for (int j = 0; j < visitedVector.size(); j++) {
        visited[visitedVector[j]] = 0;
    }

    q->clear();

    sort(nodecoll.begin(), nodecoll.end(), NodeRank::compare_level);
//	cout<<nodecoll.size()<<endl;
    int time_cnt = 0;
    for (NodeRank noderank : nodecoll) {
        int u = noderank.id;
        long start = clock();
        update_node_downward_node(k, u, q, q1, dis, dis1, visited, colors, marked_first, rev_shortcuts_arr);
        // after downhill search to fill in the objects, the state becomes stable, therefore setting marked_first[u]=0
        marked_first[u] = 0;
        time_cnt += clock() - start;
        if (time_cnt > durationLimitsMicro) {
            return 0;
        }
    }
    return 1;
}


int remove_kdnn_object_multi(int remove_node, int k_max, int k, mem_struct &mems,
                       vector<KNode> *&rev_shortcuts_arr,
                       int durationLimits, vector<KNode>* hier_local_knn_arr) {
    int durationLimitsMicro = durationLimits * 1000000;
    vector<int> visitedVector;
    DijkstraQueue *q = mems.q;
    DijkstraQueue *q1 = mems.q1;
    long *dis = mems.dist;
    long *dis1 = mems.dist1;
    int *visited = mems.visited;
    int *marked_first = mems.marked_first;
    int *colors = mems.visited1;
    dis[remove_node] = 0;
    q->push(remove_node, 0);
    visited[remove_node] = 1;
    visitedVector.push_back(remove_node);

    vector<NodeRank> nodecoll;
    while (!q->empty()) {
        int u = q->pop();
        if (u > test_n)
            continue;

        vector<KNode> &knodes = hier_local_knn_arr[u];
        vector<KNode>::iterator it;
        int is_erase = 0;
        for (it = knodes.begin(); it != knodes.end(); it++) {
            if (remove_node == it->id) {
                it = knodes.erase(it);
                is_erase = 1;
                break;
            }
        }
        // If originally, the size of the object list is less than k, that means, after object removal, it remains
        // less than k-1 objects, then all the downhill objects are already in the list. Hence, we only push in the node
        // when knodes.size()==k-1
        // set marked_first[u]=1, means that the current list is unstable, due to the object removal
        if (is_erase && knodes.size() == k - 1) {
            nodecoll.push_back(NodeRank(node_order[u], inlev_order[u], u));
            marked_first[u] = 1;
        }


        vector<KNode> *nei;
        if (mode == QUERY_FAVOR_MODE) {
            if (node_order[u] >= cut_level) {
                nei = &(level_scs_arr[cut_level][u]);
                /*
                 * if node u'level >= cut_level, then associate the node u with a level $cut_level, and employ its level
                 * shortcuts.
                 */
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
                /*
                 * if at least tradeoff_level, then "level+up"
                 */
            else
                nei = &(shortcuts_inc_arr[u]);
            /*
             * if less that tradeoff_level, then "up"
             */
        } else {
            if (node_order[u] >= cut_level) {
                continue;
            } else if (node_order[u] >= tradeoff_level)
                nei = &(shortcuts_arr[u]);
            else
                nei = &(shortcuts_inc_arr[u]);
        }


        for (KNode knode : *nei) {
            int v = knode.id;
//				int level_gap_for_v = level_gap[node_order[v]] * 2 + 1;
//				int coor_v_x = cell_loc_x[v];
//				int coor_v_y = cell_loc_y[v];
//				int gap1 = max(abs(coor_x1 - coor_v_x),
//						abs(coor_y1 - coor_v_y));
//				int gap2 = max(abs(coor_x2 - coor_v_x),
//						abs(coor_y2 - coor_v_y));
//				if (gap1 > level_gap_for_v && gap2 > level_gap_for_v)
//					continue;
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
    for (int j = 0; j < visitedVector.size(); j++) {
        visited[visitedVector[j]] = 0;
    }

    q->clear();

    sort(nodecoll.begin(), nodecoll.end(), NodeRank::compare_level);
//	cout<<nodecoll.size()<<endl;
    int time_cnt = 0;
    for (NodeRank noderank : nodecoll) {
        int u = noderank.id;
        long start = clock();
        update_node_downward_node_multi(k, u, q, q1, dis, dis1, visited, colors, marked_first,
                                        rev_shortcuts_arr, hier_local_knn_arr);
        // after downhill search to fill in the objects, the state becomes stable, therefore setting marked_first[u]=0
        marked_first[u] = 0;
        time_cnt += clock() - start;
        if (time_cnt > durationLimitsMicro) {
            return 0;
        }
    }
    return 1;
}

void remove_hier_local_knn_carpos(CarPos cp, int k_max, int k, unordered_map<int, vector<PointCar> > &car_hash,
                                  DijkstraQueue *q, long *dis, int *visited) {
    vector<int> visitedVector;
    int coor_x1 = cell_loc_x[cp.node1];
    int coor_y1 = cell_loc_y[cp.node1];
    int coor_x2 = cell_loc_x[cp.node2];
    int coor_y2 = cell_loc_y[cp.node2];
    dis[cp.node1] = cp.dis1;
    q->push(cp.node1, cp.dis1);
    visited[cp.node1] = 1;
    visitedVector.push_back(cp.node1);

    dis[cp.node2] = cp.dis2;
    q->push(cp.node2, cp.dis2);
    visited[cp.node2] = 1;
    visitedVector.push_back(cp.node2);
    vector<NodeRank> nodecoll;
    while (!q->empty()) {
        int u = q->pop();
        if (u > test_n)
            continue;

        vector<KNode> knodes = hier_local_knn_arr[u];
        vector<KNode>::iterator it;
        int is_erase = 0;
        for (it = knodes.begin(); it != knodes.end(); it++) {
            if (test_n + cp.id == it->id) {
                it = knodes.erase(it);
                is_erase = 1;
                break;
            }
        }
        hier_local_knn_arr[u] = knodes;
        if (is_erase) {
            nodecoll.push_back(NodeRank(node_order[u], inlev_order[u], u));
        }

        int need_expand = 1;
        vector<KNode> nei = shortcuts_inc_arr[u];
//		for (KNode knode : nei) {
//			int v = knode.id;
//			if (visited[v] && dis[v] + knode.dis < dis[u]) {
//				dis[u] = dis[v] + knode.dis;
//				need_expand = 0;
//				break;
//			}
//		}

        if (need_expand) {
            for (KNode knode : nei) {
                int v = knode.id;
                int level_gap_for_v = level_gap[node_order[v]] * 2 + 1;
                int coor_v_x = cell_loc_x[v];
                int coor_v_y = cell_loc_y[v];
                int gap1 = max(abs(coor_x1 - coor_v_x), abs(coor_y1 - coor_v_y));
                int gap2 = max(abs(coor_x2 - coor_v_x), abs(coor_y2 - coor_v_y));
                if (gap1 > level_gap_for_v && gap2 > level_gap_for_v)
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
    for (int j = 0; j < visitedVector.size(); j++) {
        visited[visitedVector[j]] = 0;
    }

    q->clear();
    sort(nodecoll.begin(), nodecoll.end(), NodeRank::compare);
    for (NodeRank noderank : nodecoll) {
        int u = noderank.id;
        vector<PointCar> pcs = car_hash[u];
        for (PointCar tmppc : pcs) {
            hier_local_knn_arr[u] = updateToList(hier_local_knn_arr[u], tmppc.carid + test_n, tmppc.dis, k_max);
        }
        hier_local_knn_arr[u] = update_node_downward(k_max, u, cp.id + test_n);
    }
}

void test_hierSP_correct() {
    read_road_network();
    cout << "start reading shortcuts..." << endl;
    read_scob_index();
    cout << "start testing TestShortestPath correctness..." << endl;
    mem_struct mems;
    allocate_mem(mems, test_n);

    long time1 = 0, time2 = 0;
    vector<int> level1;
    int levnum = 0;
    for (int i = 0; i < test_n; i++) {
        if (node_order[i] >= levnum) {
            level1.push_back(i);
        }
    }
    for (int i = 0; i < 5000; i++) {
        int s = get_rand_node();
        int t = get_rand_node();
        s = s % level1.size();
        t = t % level1.size();

        if (s == t)
            continue;
        s = level1[s];
        t = level1[t];
        long start1 = clock();
        int dist1 = get_distance_st(s, t);
        time1 += clock() - start1;
        long start2 = clock();
        int dist2 = TestShortestPath(s, t, node_order, shortcuts, mems,
                                     INT_MAX);
        time2 += clock() - start2;
        cout << dist1 << " " << dist2 << endl;
        if (dist1 != dist2) {
            cout << "not equal!!" << endl;
            cout << s << " " << t << endl;
            cout << dist1 << " " << dist2 << endl;

            stop_here();
        }
    }

    cout << time1 << " " << time2 << endl;
    delete_mems(mems);
}

#endif /* SRC_HIERINDEX_H_ */
