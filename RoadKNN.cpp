// AreaIntersection.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include<queue>
#include "MultiCore/RandomQueryingThread.h"
//#include "MultiCore/PartitionQueryingThread.h"
#include "data_format/MetisFormat.h"


// include order is important
#include "para/GlobalVariables.h"
#include "RoadKNNUtilities.h"
#include "dijk/NaiveKNN.h"
#include "scob/ScobIndex.h"
#include "util/DataFunctions.h"
#include "testing/TestToainFreely.h"
#include "training/ToainTraining.h"
#include "config/Config.h"
#include "training/ToainSimulation.h"
#include "training/MultiCoreTrain.h"
//#include "VTree-master/vtree_build.h"


using namespace std;

extern vector<KNode> **shortcuts_arr_local;
extern vector<int> **cells;
extern int test_n;
extern int cell_num;
//3000000 for NY, FLA
//7000000 for W
extern int *regional_landmarks_all;

//deprecated
int get_grid_dis(int s, int t) {
    int s_cell_x = cell_loc_x[s];
    int t_cell_x = cell_loc_x[t];
    int s_cell_y = cell_loc_y[s];
    int t_cell_y = cell_loc_y[t];
    return max(abs(s_cell_x - t_cell_x), abs(s_cell_y - t_cell_y));
}


/**
 * We need to initialize the vectors, arrays, etc.
 */

// deprecated
void refresh_cars() {
    read_road_network();

    // build a grid environment for estimate the nearest node of the object
    int grid_resolution_for_estimate_car_positions = 500;
    cell_num = grid_resolution_for_estimate_car_positions;
    cells = new vector<int> *[cell_num];
    circle_cells = new vector<int> *[cell_num];
    for (auto i = 0; i < cell_num; ++i) {
        cells[i] = new vector<int>[cell_num];
        circle_cells[i] = new vector<int>[cell_num];
    }
    cout << "start build cell..." << endl;
    build_cells();

    cout << "start get car position..." << endl;
    if (data_space.compare("0") == 1)
        get_car_position_beijing();
    if (data_space.compare("1") == 1) {
        //	get_car_position_ny();
        get_car_position_random();
    }

    // dispose
    for (auto i = 0; i < cell_num; ++i) {
        delete[] cells[i];
        delete[] circle_cells[i];
    }
    delete[] cells;
    delete[] circle_cells;
}


void cal_rt(const char *path, const char *path_out, int test_mode) {
    // test_mode=1 uber
    // test_mode=0 poke
    FILE *file = fopen((input_parameters.input_data_dir + path).c_str(), "r");
    if (file == NULL) {
        cout << "file of " << path << " not found, in cal_rt" << endl;
        stop_here();
    }
    char title1[20], title2[20], title3[20], title4[20], title5[20], title6[20];
    fscanf(file, "%s %s %s %s %s %s", title1, title2, title3, title4, title5, title6);
    cout << title3 << endl;
    FILE *fout = fopen((input_parameters.output_data_dir + path_out).c_str(), "w");
    fprintf(fout, "%s \t %s \t %s \t %s \t %s \t %s \t R_t \t T \t RT\n", title1, title2, title3, title4, title5,
            title6);
    double T_q, V_q, T_u, V_u, ratio;
    int k;
    double microsec_per_sec = 1000000.0;
    if (test_mode == 1) {
        while (fscanf(file, "%d %lf %lf %lf %lf %lf", &k, &ratio, &T_q, &V_q, &T_u, &V_u) == 6) {
            cout << k << "\t" << ratio << "\t" << T_q << "\t" << V_q << "\t" << T_u << "\t" << V_u << endl;
            for (double T = 1.0; T <= 16.0; T *= 2) {
                for (double T_r = 100.0; T_r < 10000.0; T_r *= 2) {
                    double rt = analytical_query_first(T_q / microsec_per_sec, T_u / microsec_per_sec,
                                                       T_r / microsec_per_sec, test_n * ratio, T,
                                                       V_q / microsec_per_sec / microsec_per_sec);

                    fprintf(fout, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", k, ratio, T_q, V_q, T_u, V_u, T_r, T, rt);

                }
            }

            //		cal_poke_RT(double T_q, double T_u, double T_r, double lambda_u, double V_q, double V_u);
        }
    } else {
        while (fscanf(file, "%d %lf %lf %lf %lf %lf", &k, &ratio, &T_q, &V_q, &T_u, &V_u) == 6) {
            cout << k << "\t" << ratio << "\t" << T_q << "\t" << V_q << "\t" << T_u << "\t" << V_u << endl;
            for (double T = 1.0; T <= 16.0; T *= 2) {
                for (double T_r = 100.0; T_r < 10000.0; T_r *= 2) {
                    double rt = analytical_fifo(T_q / microsec_per_sec, T_u / microsec_per_sec,
                                                T_r / microsec_per_sec, test_n * ratio / T,
                                                V_q / microsec_per_sec / microsec_per_sec,
                                                V_u / microsec_per_sec / microsec_per_sec);
                    fprintf(fout, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", k, ratio, T_q, V_q, T_u, V_u, T_r, T, rt);

                }
            }

        }

    }
    fclose(file);
    fclose(fout);

}


int pair_compare_second(const std::pair<int, int> &left, const std::pair<int, int> &right) {
    return left.second > right.second;
}

vector<std::pair<int, int> > gen_test_from_taxi_data() {
    read_road_network();
    vector<int> start_nodes;
    vector<int> end_nodes;
    int *frequency = new int[test_n + 1];
    memset(frequency, 0, sizeof(int) * (test_n + 1));
    char buffer1[100], buffer2[100], buffer3[100], buffer4[100];
    int d1 = 1, d2 = 1;
    char c1 = 'c', c2 = 'c';
    FILE *fInput = fopen((input_parameters.input_data_dir + "yellow_tripdata_2015-01_full.csv").c_str(), "r");
    if (fInput == NULL) {
        cout << "fInput==NULL" << endl;
        char c;
        cin >> c;
    }
    double x1 = 1, y1 = 1, x2 = 1, y2 = 1, n1 = 1, n2 = 1;
    int pcnt = 1;
    double nw_dist = 1;
    fscanf(fInput, "%*s");
    int index = 0;
    int cnt = 0;
//    while (fscanf(fInput, "%s %s %13s%lf,%lf,%lf,%lf,%d,%lf,%s", buffer1,
//                  buffer2, buffer3, &x1, &y1, &x2, &y2, &pcnt, &nw_dist, buffer4)
//           == 10) {
    while (fscanf(fInput, "%s %s %8s,%d,%lf,%lf,%lf,%d,%c,%lf,%lf,%s", buffer1, buffer2, buffer3, &pcnt, &nw_dist, &x1,
                  &y1, &d2, &c1, &x2, &y2, buffer4) == 12) {
        cout << cnt << endl;
        if (cnt > 1000000)
            break;
//        if (nw_dist < 5)
//            continue;
        cnt++;
        int node1 = getNearestNode(x1, y1);
        int node2 = getNearestNode(x2, y2);
        //		cout << node1 << " " << node2 << endl;
        frequency[node1]++;
        frequency[node2]++;
//        start_nodes.push_back(node1);
//        end_nodes.push_back(node2);
    }

    fclose(fInput);
    vector<std::pair<int, int> > node_ranks;
    for (int i = 0; i < test_n; i++) {
        node_ranks.push_back(make_pair(i, frequency[i]));
    }
    std::sort(node_ranks.begin(), node_ranks.end(), pair_compare_second);
    FILE *node_frequency_ranks_file = fopen((input_parameters.input_data_dir + "node_freq_ranks.txt").c_str(), "w");
    for (pair<int, int> p: node_ranks) {
        fprintf(node_frequency_ranks_file, "%d %d\n", p.first, p.second);
    }
    fclose(node_frequency_ranks_file);
    return node_ranks;
}

vector<pair<int, int> > read_taxi_ranks() {
    FILE *file = fopen((input_parameters.input_data_dir + "node_freq_ranks.txt").c_str(), "r");
    cout << (input_parameters.input_data_dir + "node_freq_ranks.txt").c_str() << endl;
//    char c; cin>>c;
    test_valid_file(file, (input_parameters.input_data_dir + "node_freq_ranks.txt").c_str(), "read_taxi_ranks");
    int node, rank;
    while (fscanf(file, "%d %d\n", &node, &rank) == 2) {
        node_ranks.push_back(make_pair(node, rank));
    }
    fclose(file);
    return node_ranks;
}

void generate_order_from_betweenness() {
    FILE *in_file = fopen((input_parameters.input_data_dir + "appr_between_order.txt").c_str(), "r");
    test_valid_file(in_file, (input_parameters.input_data_dir + "appr_between_order.txt"),
                    "generate_order_from_betweenness");
    string file_name1 = data_space;
    file_name1 = file_name1.append("between_node_order.txt");
    file_name1 = input_parameters.input_data_dir + file_name1;
    FILE *nodeorder_file;
    nodeorder_file = fopen(file_name1.c_str(), "w");
    vector<int> ids;
    vector<double> rank_values;
    int id;
    double rv;
    int h = 7;
    double max_val = -1.0;
    double min_val = 100000000.0;
    while (fscanf(in_file, "%d %lf\n", &id, &rv) == 2) {
        ids.push_back(id);
        //if(rv<0) rv=0;
        rank_values.push_back(rv);
        if (rv > max_val) max_val = rv;
        if (rv < min_val) min_val = rv;
    }
    double gap = (max_val - min_val) / h + 1;
    cout << max_val << endl;
    cout << min_val << endl;
    cout << "gap: " << gap << endl;
    for (int i = 0; i < ids.size(); i++) {
        fprintf(nodeorder_file, "%d %d\n", ids[i] + 1, (int) ((rank_values[i] - min_val) / gap));
    }


    fclose(nodeorder_file);
    fclose(in_file);
}

void generate_order_from_coverage_number() {

    FILE *in_file = fopen((input_parameters.input_data_dir + "adaptive_coverage_output_order.txt").c_str(), "r");
    test_valid_file(in_file, (input_parameters.input_data_dir + "adaptive_coverage_output_order.txt"),
                    "generate_order_from_coverage");
    string file_name1 = data_space;
    file_name1 = file_name1.append("adap_coverage_node_order.txt");
    file_name1 = input_parameters.input_data_dir + file_name1;
    FILE *nodeorder_file;
    nodeorder_file = fopen(file_name1.c_str(), "w");
    vector<int> ids;
    vector<int> rank_values;
    int id, rv;
    int h = 10;
    int max_val = -1;
    while (fscanf(in_file, "%d %d\n", &id, &rv) == 2) {
        ids.push_back(id);
        rank_values.push_back(rv);
        if (rv > max_val) max_val = rv;
    }


    vector<std::pair<int, int> > node_ranks;
    for (int i = 0; i < ids.size(); i++) {
        node_ranks.push_back(make_pair(ids[i] + 1, rank_values[i]));
    }
    std::sort(node_ranks.begin(), node_ranks.end(), pair_compare_second);
    int gap_size = ids.size() / h + 1;
    for (int i = 0; i < node_ranks.size(); i++) {
        fprintf(nodeorder_file, "%d %d\n", node_ranks[i].first, int(i / gap_size));
    }

    fclose(nodeorder_file);
    fclose(in_file);
    cout << "coverage generaged!" << endl;
}


void generate_order_from_betweenness_number() {

    FILE *in_file = fopen((input_parameters.input_data_dir + "adaptive_between_output_order.txt").c_str(), "r");
    test_valid_file(in_file, (input_parameters.input_data_dir + "adaptive_between_output_order.txt"),
                    "generate_order_from_betweenness");
    string file_name1 = data_space;
    file_name1 = file_name1.append("adap_between_node_order.txt");
    file_name1 = input_parameters.input_data_dir + file_name1;
    FILE *nodeorder_file;
    nodeorder_file = fopen(file_name1.c_str(), "w");
    vector<int> ids;
    vector<double> rank_values;
    int id;
    double rv;
    int h = 10;
    while (fscanf(in_file, "%d %lf\n", &id, &rv) == 2) {
        ids.push_back(id);
        rank_values.push_back(rv);

    }


    vector<std::pair<int, int> > node_ranks;
    for (int i = 0; i < ids.size(); i++) {
        node_ranks.push_back(make_pair(ids[i] + 1, rank_values[i]));
    }
    std::sort(node_ranks.begin(), node_ranks.end(), pair_compare_second);
    int gap_size = ids.size() / h + 1;
    for (int i = 0; i < node_ranks.size(); i++) {
        fprintf(nodeorder_file, "%d %d\n", node_ranks[i].first, i / gap_size);
    }

    fclose(nodeorder_file);
    fclose(in_file);
    cout << "betweenness gernated!" << endl;
}

void generate_order_from_coverage() {

    FILE *in_file = fopen((input_parameters.input_data_dir + "appr_coverage_order.txt").c_str(), "r");
    test_valid_file(in_file, (input_parameters.input_data_dir + "appr_coverage_order.txt"),
                    "generate_order_from_coverage");
    string file_name1 = data_space;
    file_name1 = file_name1.append("coverage_node_order.txt");
    file_name1 = input_parameters.input_data_dir + file_name1;
    FILE *nodeorder_file;
    nodeorder_file = fopen(file_name1.c_str(), "w");
    vector<int> ids;
    vector<int> rank_values;
    int id, rv;
    int h = 7;
    int max_val = -1;
    while (fscanf(in_file, "%d %d\n", &id, &rv) == 2) {
        ids.push_back(id);
        rank_values.push_back(rv);
        if (rv > max_val) max_val = rv;
    }
    int gap = max_val / h + 1;
    cout << max_val << endl;
    for (int i = 0; i < ids.size(); i++) {
        fprintf(nodeorder_file, "%d %d\n", ids[i] + 1, rank_values[i] / gap);
    }


    fclose(nodeorder_file);
    fclose(in_file);
}


int main(int argc, char *argv[]) {
    // demo code
    // Example procedure:
    // set

//    cout << argc << endl;


    if (argc % 2 == 0) {
        cout << "args not correct!" << endl;
        stop_here();
    }
    for (int i = 1; i < argc; i += 2) {
        std::string str(argv[i]);
        paras[str] = argv[i + 1];
    }
    multiTestPara = anaylizeParameters();
    srand(time(NULL));
    set_test_para();
    read_road_network();

    if (paras.find("-multicore") != paras.end()) {
        struct timeval start, end;
        srand(time(NULL));

        double query_rate = multiTestPara.query_rate/multiTestPara.layer;
        double insert_rate = multiTestPara.insert_rate;
        double delete_rate = multiTestPara.delete_rate;
        multiTestPara.num_total_threads/=multiTestPara.layer;
        double alpha = 1.0;
        int k = 10;
        double fail_p = multiTestPara.pfail;
        can_estimate = 1;


        cout << "num of threads: " << multiTestPara.num_threads_update << endl;
        string loadfile;
        if (multiTestPara.method_name.compare("vtree") == 0) {
            car_id_insert_cnt =  new int[MAX_CORES];

            memset(car_id_insert_cnt, 0, sizeof(int)*MAX_CORES);

            loadfile = input_parameters.input_data_dir + network_name + ".vtree";
            if(network_name.compare("BJ-old")==0){
                loadfile = input_parameters.input_data_dir  + "BJ.vtree";

            }
            cout << "loadfile: " << loadfile << endl;
            init();
            if(multiTestPara.auto_config)
                ChooseVtreePRConfiguration(fail_p, alpha, k, multiTestPara);

        }
        int configurationId = 9;
        cout<<"method name: "<<multiTestPara.method_name<<endl;

        if (multiTestPara.method_name.compare("toain") == 0) {
            init_memory();
            cout<<"network_name: "<<network_name<<endl;
            cout<<"toain type: "<<multiTestPara.toain_type<<endl;
            configurationId = chooseSingleTOAINConfiguration(fail_p, alpha, k,
                                                             multiTestPara);
            cout<<"configurationID: "<<configurationId<<endl;
//            multiTestPara = anaylizeParameters();
            multiTestPara.num_total_threads/=multiTestPara.layer;
            multiTestPara.num_threads_query=get_num_threads_query(multiTestPara.num_threads_update, multiTestPara.is_single_aggregate,
                                                                  multiTestPara.num_total_threads);
            if(multiTestPara.auto_config) {
                ChooseTOAINPRConfiguration(fail_p, alpha, k,
                                           configurationId, multiTestPara);
            }
            else {
                prepare_environment_with_scob_conf(configurationId);
            }
	        multiTestPara.num_threads_query = multiTestPara.query_thread;
            multiTestPara.num_threads_update = multiTestPara.update_thread;
            cout << "finish choosing configuration!" << endl;
            cout << "multiTestPara.num_threads_query: " << multiTestPara.num_threads_query << endl;
            cout << "multiTestPara.num_threads_update: " << multiTestPara.num_threads_update << endl;
            cout << "start online running" << endl;
        }
        if (multiTestPara.method_name.compare("dijk") == 0) {
            if(multiTestPara.auto_config) {
                ChooseDIJKPRConfiguration(fail_p, alpha, k, multiTestPara);
            }
            else{
                multiTestPara.num_threads_query = multiTestPara.query_thread;
                multiTestPara.num_threads_update = multiTestPara.update_thread;
                cout << "finish choosing configuration!" << endl;
                cout << "multiTestPara.num_threads_query: " << multiTestPara.num_threads_query << endl;
                cout << "multiTestPara.num_threads_update: " << multiTestPara.num_threads_update << endl;
                cout << "start online running" << endl;
            }
        }

        int test_num = 1;
        double response_time_all = 0.0;
        double query_finish_rate_all = 0.0;
        double query_process_time_all = 0.0;
        double update_process_time_all = 0.0;
        double schedule_cost_all = 0.0;
        for(int i=0;i<test_num;i++) {
            if (multiTestPara.method_name.compare("vtree") == 0) {
                vtrees.clear();
                cout << "total index copies: " << multiTestPara.num_threads_query * multiTestPara.num_threads_update
                     << endl;
                for (int i = 0; i < multiTestPara.num_threads_query * multiTestPara.num_threads_update; i++) {
                    cout << "loading " << i << endl;
                    load_binary(loadfile, i);
                    vtrees.push_back(&(tree[i]));
                }
            }
            double object_ratio = 0.001;
            if (multiTestPara.num_threads_update != -1) {
                if (multiTestPara.parmethod.compare("rand") == 0) {
                    overload_flag = 0;

                    is_simulation = 1;
                    RandomThreadPool_Control *tp = new RandomThreadPool_Control(0, 0, test_n, multiTestPara.num_total_threads, alpha, k, fail_p,test_n,query_rate, insert_rate,
                                                                                delete_rate, multiTestPara.test_simulation_time,multiTestPara.query_cost,
                                                                                multiTestPara.insert_cost,multiTestPara.delete_cost,configurationId);
                    tp->run();


//            RandomThreadPool *tp2 = new RandomThreadPool(0, test_n/2, test_n, multiTestPara.num_threads_query,
//                                                         multiTestPara.num_threads_update, alpha, k, fail_p, test_n,
//                                                         query_rate, insert_rate,
//                                                         delete_rate, multiTestPara.test_simulation_time);

                    }

                } else {

                }





            }
            number_of_queries = 1;
            total_update_response_time = 0.0;
            number_of_updates = 1;
            total_update_process_time = 0.0;
            total_query_process_time = 0.0;
            number_of_query_processings = 1;
        }
//        std::ofstream outfile;
//
//        outfile.open(input_parameters.output_data_dir + "outfile" + (multiTestPara.suffix), std::ios_base::app);
//
//        outfile << endl
//                << network_name<<" "
//                << "init: "<<multiTestPara.init_objects<<" "
//                << multiTestPara.method_name << " config simulation time: "
//                << multiTestPara.config_simulation_time << " test simulate time: "
//                << multiTestPara.test_simulation_time << " configure: "
//                << configurationId << " threshold: " << multiTestPara.is_thresholded << " fail_p: " << fail_p
//                << " "
//                << query_rate << " " << insert_rate
//                << " " << delete_rate << " " << multiTestPara.method_name << " singleAggregate: "
//                << multiTestPara.is_single_aggregate << " "
//                << multiTestPara.num_threads_update
//                << " " << multiTestPara.num_threads_query << " "
//                <<"average response time: "<<response_time_all/test_num
//                <<" average query time: "<<query_process_time_all / test_num
//                <<" average update time: "<<update_process_time_all / test_num
//                <<" average schedule cost: "<<schedule_cost_all / test_num
//                <<" average query finish: "<<query_finish_rate_all/test_num
//                <<" layer: "<<multiTestPara.layer
//                << endl;
//
//
//        outfile.close();

        if (multiTestPara.method_name.compare("toain") == 0) {
            delete_memory();
        }



        cout << endl << endl << endl << endl << endl;




    if (test_parameters.demo) {

        init_memory();
        read_scob_index();
        cout << "the number of nodes: " << test_n << endl;
        level_num = get_level_num_with_dataspace();

        vector<ScobConfiguration> configurations = ScobConfiguration::generate_all_configurations(
                level_num);
        int test_case = 0;
        tradeoff_level = configurations[test_case].tradeoff_level;
        cut_level = configurations[test_case].cut_level;
        mode = configurations[test_case].mode;
        cout << "tradeoff_level: " << tradeoff_level << endl;
        cout << "cut_level: " << cut_level << endl;
        cout << "mode: " << mode << endl;
        int k = 10;
        int k_max = 50;// default value for k_max
        mem_struct mems;
        allocate_mem(mems, test_n + test_n + 1);
        int *car_nodes = new int[test_n];
        memset(car_nodes, 0, test_n * sizeof(int));
        objects_vector.clear();
        // build a scob with 1000 objects
        for (int i = 0; i < 1000; i++) {
            int object_node = get_rand_node();
            int is_init = 1;
            // insert an object
            update_kDNN_object(is_init, object_node, k, 0, mems.dist, mems.visited, mems.q);
            car_nodes[object_node] = 1;
            objects_vector.push_back(object_node);
        }
        // perform a query
        int query_node = 10000;
        int removeNum = 5;
        // perform a Dijkstra
        vector<KNode> kNNs = naiveKNN(k, objects_vector, query_node, mems.dist, mems.visited, mems.q, car_nodes);
        cout << "kNNs: " << endl;
        for (KNode knode: kNNs) {
            cout << knode.id << " " << knode.dis << endl;
        }
        for (int i = 0; i < removeNum; i++) {
            car_nodes[kNNs[i].id] = 0;
        }
        kNNs = naiveKNN(k, objects_vector, query_node, mems.dist, mems.visited, mems.q, car_nodes);
        cout << "after removals kNNs: " << endl;
        for (KNode knode: kNNs) {
            cout << knode.id << " " << knode.dis << endl;
        }
        delete[] car_nodes;

        kNNs = scob_KNN_query(k, hier_local_knn_arr, query_node, mems.dist, mems.visited, mems.q);

        cout << "kNNs: " << endl;
        for (KNode knode: kNNs) {
            cout << knode.id << " " << knode.dis << endl;
        }



        // remove objects, with performance guarantee as 10 seconds
        vector<int> removedObjectVec;

        for (int i = 0; i < removeNum; i++) {
            int removeObject = kNNs[i].id;
            removedObjectVec.push_back(removeObject);
            cout << "remove object: " << removeObject << endl;
            int success_remove = remove_kdnn_object(removeObject, k_max, k, mems, rev_shortcuts_arr, 10);

            cout << "after removal, kNNs: " << endl;
            if (success_remove) cout << "removal successful" << endl;
            else cout << "removal unsuccessful" << endl;
        }
        kNNs = scob_KNN_query(k, hier_local_knn_arr, query_node, mems.dist, mems.visited, mems.q);
        for (KNode knode: kNNs) {
            cout << knode.id << " " << knode.dis << endl;
        }
        // perform inserts
        for (int i = 0; i < removeNum; i++) {
            int removeObject = removedObjectVec.at(i);
            int is_init = 0;
            update_kDNN_object(is_init, removeObject, k_max, k, mems.dist, mems.visited, mems.q);
        }
        cout << "after insert" << endl;
        cout << "kNNs: " << endl;
        kNNs = scob_KNN_query(k, hier_local_knn_arr, query_node, mems.dist, mems.visited, mems.q);
        for (KNode knode: kNNs) {
            cout << knode.id << " " << knode.dis << endl;
        }


        delete_mems(mems);
        delete_memory();


    }
//    stop_here();
    return 0;

}

