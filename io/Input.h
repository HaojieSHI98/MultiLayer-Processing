//
// Created by lostrong on 17/7/10.
//

#ifndef SOB_INPUT_H
#define SOB_INPUT_H

#include <iostream>
#include "../para/GlobalVariables.h"
#include "../util/GeneralUtility.h"
using namespace std;
FILE* fedge, *fnode;
extern int eNum;
extern int eNext[2 * MAGIC_NUM + 2], eDis[2 * MAGIC_NUM + 1], eNode[2 * MAGIC_NUM + 1];
void get_input(string filenode, string fileedge) {
    cout << "reading files: " << filenode << " " << fileedge << endl;

    x = new int[MAGIC_NUM];
    y = new int[MAGIC_NUM];

    cell_loc_x = new int[MAGIC_NUM];
    cell_loc_y = new int[MAGIC_NUM];
    memset(cell_loc_x, 0, sizeof(int) * MAGIC_NUM);
    memset(cell_loc_y, 0, sizeof(int) * MAGIC_NUM);

    fnode = fopen((filenode).c_str(), "r");
    fedge = fopen((fileedge).c_str(), "r");

    test_valid_file(fnode, (filenode), "get_input");
    test_valid_file(fnode, (fileedge), "get_input");


    cout << "scanning node files ..." << endl;
    int i = 0;
    int s, t, d;
    double doubleD, doubleT;
if (network_map[input_parameters.graph_node_file].compare("BJ-old")==0){
    while (fscanf(fnode, "%d,%lf,%lf\n", &s, &doubleD, &doubleT) == 3) {
        x[i] = (int) (doubleT * 1000000);
        y[i] = (int) (doubleD * 1000000);
        ++i;
    }
}
    else {
    while (fscanf(fnode, "v %d%d%d\n", &s, &t, &d) == 3) {
        x[i] = t;
        y[i] = d;
        if(x_min> x[i]) x_min = x[i];
        if(x_max< x[i]) x_max = x[i];
        if(y_min > y[i]) y_min = y[i];
        if(y_max < y[i]) y_max = y[i];
        ++i;
    }
}
    fclose(fnode);
    test_n = i;

    regional_landmarks_all = new int[test_n];
    cell_loc = new int[test_n + 1];
    for (i = 0; i < NUM_CELL_IN_NEIGHBOR; ++i) {
        circle_cell_loc[i] = new int[test_n + 1];
    }
    vStart = new int[test_n + 1];

    memset(vStart, 0, sizeof(int) * (test_n + 1));

    i = 0;
    int test_data_cnt = 0;
    int test_prev_d = -1;
    int sum_edge = 0;
    int pivot_edge = 100000;
    int pivot_cnt = 0;
    cout << "scanning edge files" << endl;
    if (network_map[input_parameters.graph_node_file].compare("BJ-old")==0){
        int edgeId;
        double doubleEdgeD;
        while (fscanf(fedge, "%d,%d,%d,%lf\n", &edgeId, &s, &t, &doubleEdgeD) == 4) {
            ++i;
            valid_node_hash[s] = 1;
            valid_node_hash[t] = 1;
            eNext[eNum] = vStart[s];
            eDis[eNum] = (int) (doubleEdgeD * 1000000);
            eNode[eNum] = t;
            vStart[s] = eNum++;
            sum_edge += get_eu_dist(s, t);
            pivot_cnt += sum_edge / pivot_edge;
            sum_edge %= pivot_edge;
            if ((int) (doubleEdgeD * 1000000) != test_prev_d) {
                if (test_data_cnt % 2 != 0) {
                    cout << "data error: not undirected!" << endl;
                    //			cout << s << " " << t << " " << (int)(doubleEdgeD*1000000) << endl;
                    stop_here();
                }

                test_data_cnt++;
                test_prev_d = (int) (doubleEdgeD * 1000000);
            } else {
                test_data_cnt--;
            }
        }
    }
    else {
        while (fscanf(fedge, "a %d%d%d\n", &s, &t, &d) == 3) {
            ++i;
            valid_node_hash[s] = 1;
            valid_node_hash[t] = 1;
            eNext[eNum] = vStart[s];
            eDis[eNum] = d;
            eNode[eNum] = t;
            vStart[s] = eNum++;
            sum_edge += get_eu_dist(s, t);
            pivot_cnt += sum_edge / pivot_edge;
            sum_edge %= pivot_edge;
            if (d != test_prev_d) {
                if (test_data_cnt % 2 != 0) {
                    cout << "data error: not undirected!" << endl;
                    cout << s << " " << t << " " << d << endl;
                    stop_here();
                }

                test_data_cnt++;
                test_prev_d = d;
            } else {
                test_data_cnt--;
            }
        }
    }
    test_m = i;
    fclose(fedge);

    cout << "avg edge weight:" << pivot_cnt * (pivot_edge * 1.0 / test_m) << endl;
    cout << "#nodes: " << test_n << endl;
    cout << "#edges: " << test_m << endl;
}

// deprecated
// for shenzhou dataset
void get_input_shenzhou(string filenode, string fileedge) {
    data_space = data_map[filenode];
    network_name = network_map[filenode];
    cout << "reading files: " << filenode << " " << fileedge << endl;

    x = new int[25000000];
    y = new int[25000000];

    cell_loc_x = new int[MAGIC_NUM];
    cell_loc_y = new int[MAGIC_NUM];
    memset(cell_loc_x, 0, sizeof(int) * MAGIC_NUM);
    memset(cell_loc_y, 0, sizeof(int) * MAGIC_NUM);

    fnode = fopen((filenode).c_str(), "r");
    fedge = fopen((fileedge).c_str(), "r");
    int s, t;
    double doubleD, doubleT;
    test_valid_file(fnode, (filenode), "get_input");
    test_valid_file(fnode, (fileedge), "get_input");


    cout << "scanning node files ..." << endl;
    int i = 0;
    while (fscanf(fnode, "%d,%lf,%lf\n", &s, &doubleD, &doubleT) == 3) {
        x[i] = (int) (doubleT * 1000000);
        y[i] = (int) (doubleD * 1000000);
        ++i;
    }
    fclose(fnode);
    test_n = i;

    regional_landmarks_all = new int[test_n];
    cell_loc = new int[test_n + 1];
    for (i = 0; i < NUM_CELL_IN_NEIGHBOR; ++i) {
        circle_cell_loc[i] = new int[test_n + 1];
    }
    vStart = new int[test_n + 1];

    memset(vStart, 0, sizeof(int) * (test_n + 1));

    i = 0;
    int test_data_cnt = 0;
    int test_prev_d = -1;
    int sum_edge = 0;
    int pivot_edge = 100000;
    int pivot_cnt = 0;
    cout << "scanning edge files" << endl;
    int edgeId;
    double doubleEdgeD;
    while (fscanf(fedge, "%d,%d,%d,%lf\n", &edgeId, &s, &t, &doubleEdgeD) == 4) {
        ++i;
        valid_node_hash[s] = 1;
        valid_node_hash[t] = 1;
        eNext[eNum] = vStart[s];
        eDis[eNum] = (int) (doubleEdgeD * 1000000);
        eNode[eNum] = t;
        vStart[s] = eNum++;
        sum_edge += get_eu_dist(s, t);
        pivot_cnt += sum_edge / pivot_edge;
        sum_edge %= pivot_edge;
        if ((int) (doubleEdgeD * 1000000) != test_prev_d) {
            if (test_data_cnt % 2 != 0) {
                cout << "data error: not undirected!" << endl;
                //			cout << s << " " << t << " " << (int)(doubleEdgeD*1000000) << endl;
               stop_here();
            }

            test_data_cnt++;
            test_prev_d = (int) (doubleEdgeD * 1000000);
        } else {
            test_data_cnt--;
        }
    }
    test_m = i;
    fclose(fedge);
    cout << "avg edge:" << pivot_cnt * (pivot_edge * 1.0 / test_m) << endl;
    cout << "test_n: " << test_n << " test_m: " << test_m << endl;

}

vector<int> read_network_parts(int parts, int num_nodes){
    vector<int> part_maps;
    if(parts>1) {
        string in_path = network_name + "_";
        in_path = in_path.append("metis_graph.txt.part.");
        stringstream ss;
        ss << parts;

        in_path = in_path.append(ss.str());
        in_path = input_parameters.input_data_dir + in_path;
        cout << in_path << endl;
        FILE *in_file = fopen(in_path.c_str(), "r");
        if (in_file == NULL) {
            cout << "in_file==NULL" << endl;
            stop_here();
        }

        int part_id;
//    int id_cnt=1;
        part_maps.push_back(0);
        while (fscanf(in_file, "%d\n", &part_id) == 1) {
            part_maps.push_back(part_id);
//        id_cnt++;
        }
    }
    else{
        for(int i=0;i<num_nodes;i++) part_maps.push_back(0);
    }
    return part_maps;

}
vector<int> get_network_parts_by_coor(int parts){
    vector<int> part_maps;
    part_maps.push_back(0);
    int x_mid = (x_max+x_min)/2;
    int y_mid = (y_max+y_min)/2;
    if(parts==4) {
        for (int i = 1; i <= test_n; i++) {
            if (x[i - 1] < x_mid && y[i - 1] < y_mid) {
                part_maps.push_back(0);
            }
            if (x[i - 1] < x_mid && y[i - 1] >= y_mid) {
                part_maps.push_back(1);
            }
            if (x[i - 1] >= x_mid && y[i - 1] < y_mid) {
                part_maps.push_back(2);
            }
            if (x[i - 1] >= x_mid && y[i - 1] >= y_mid) {
                part_maps.push_back(3);
            }
        }
    }
    if(parts==2){
        int cnt0=0, cnt1=0;
        for (int i = 1; i <= test_n; i++) {
            if (x[i - 1] < x_mid ) {
                part_maps.push_back(0);
                cnt0++;
            }

            if (x[i - 1] >= x_mid) {
                part_maps.push_back(1);
                cnt1++;
            }
        }
        cout<<cnt0<<" "<<cnt1<<endl;
//        stop_here();

    }
    if(parts==1){
        for (int i = 1; i <= test_n; i++) {
            part_maps.push_back(0);
        }

    }
    return part_maps;

}
void read_road_network() {
    if (has_read_road_network == 0) {
        cout << "Read Roadnetworks begin ..." << endl;

        get_input(input_parameters.input_data_dir+input_parameters.graph_node_file,
                  input_parameters.input_data_dir+input_parameters.graph_edge_file);
        cout << "Read Roadnetworks end ..." << endl;
        has_read_road_network = 1;
    }

}

#endif //SOB_INPUT_H
