//
// Created by lostrong on 17/3/9.
//

#ifndef SOB_TESTING_H
#define SOB_TESTING_H

#include "../util/TimeStat.h"
#include "../para/GlobalVariables.h"
#include <vector>
#include <sstream>
#include <string.h>
#include "../scob/ScobQuery.h"
#include "../util/TimeStat.h"
#include "TestToainFreely.h"
#include "../util/DataFunctions.h"
#include "../RoadKNNUtilities.h"
#include "../queue/GeneralQueueSimulator.h"


double analytical_fifo(double T_q, double T_u, double T_r, double lambda_u, double V_q, double V_u) {
//	M/M/1 return (T_r - T_q - lambda_u * T_u * T_u - lambda_u * T_u * (T_r - T_q)) / (T_q * T_r);
//	M/G/1
    if (T_r < T_q)
        return 0;
    long double x3 =
            (2 * (T_r - T_q) * (1 - lambda_u * T_u) - lambda_u * (V_u + T_u * T_u)) / (V_q + 2 * T_r * T_q - T_q * T_q);
    long double x2 = (1 - lambda_u * T_u) / T_q;
    double c1 = x3 < x2 ? x3 : x2;
    if (c1 < 0) return 0;
    return c1;
}

double analytical_query_first(double T_q, double T_u, double T_r, double tau, double T, double V_q) {
//	M/M/1 return min(1.0 / T_q - 1.0 / T_r, 1.0 / T_q - tau * T_u / (T * T_q));
//	M/G/1
    if (T_r < T_q)
        return 0;
    double a = 2 * (T_r - T_q) / (V_q + 2 * T_r * T_q - T_q * T_q);
    double lambda_u = tau / T;
    double b = (1 - lambda_u * T_u) / T_q;
    double c = a < b ? a : b;
    if (c < 0) return 0;
    return c;
}

void analytical_fifo_all_parameters(vector<times> &cgraph_times, vector<times> &dijk_times, double obj_num, FILE *file,
                                    vector<int> &accumulated_comb_ids) {
    // Calculate Pokemon-Go's RT (Response time guaranteed Throughput) for each algorithm
    // T_u of Dijkstra's algorithm is 0, T_q of Dijkstra's algorithm is dijk_query_times[i]
    // T_u of CGraph's algorithm is
    cout << "Pokemon-Go's model:" << endl;

    for (int T_r = 100; T_r <= 10000; T_r *= 2) {
        cout << "T_r is " << T_r << endl;
        for (double T = 1.0; T <= 16.0; T *= 2) {
            vector<double> Dijkstra_Thr, CGraph_Thr;
            double microsec_per_sec = 1000000.0;
            for (int alg = 0; alg < dijk_times.size(); alg++) {
                Dijkstra_Thr.push_back(
                        analytical_fifo(dijk_times[alg].query / microsec_per_sec, 0.0, T_r / microsec_per_sec,
                                        obj_num / T,
                                        dijk_times[alg].query_variance / microsec_per_sec / microsec_per_sec, 0.0));
                CGraph_Thr.push_back(
                        analytical_fifo(cgraph_times[alg].query / microsec_per_sec,
                                        cgraph_times[alg].update / microsec_per_sec, T_r / microsec_per_sec,
                                        obj_num / T,
                                        cgraph_times[alg].query_variance / microsec_per_sec / microsec_per_sec,
                                        cgraph_times[alg].update_variance / microsec_per_sec / microsec_per_sec));

            }
            double max_thr = -1;
            for (double thr : Dijkstra_Thr) {
                if (max_thr < thr) {
                    max_thr = thr;
                }
            }
            cout << "max throughput of Dijkstra's algorithm: " << max_thr << endl;

            double max_thr_ours = -1;
            int comb_id = -1;
            int ind = 0;
            for (double thr : CGraph_Thr) {
                if (max_thr_ours < thr) {
                    max_thr_ours = thr;
                    comb_id = ind;
                }
                ind++;
            }
            cout << "max throughput of CGraph's algorithm: " << max_thr_ours << endl;
            cout << "combination preferred is combination " << comb_id << endl;
            accumulated_comb_ids.push_back(comb_id);
            cout << endl;
            fprintf(file, "%f \t %f \t %d \t %f \t %f \t %d\n", obj_num / test_n, T, T_r, max_thr_ours, max_thr,
                    comb_id);
        }
    }
}

void analytical_query_first_all_parameters(vector<times> &cgraph_times, vector<times> &dijk_times, double obj_num,
                                           FILE *file,
                                           vector<int> &accumulated_comb_ids) {
    // Calculate Uber's RT (Response time guaranteed Throughput) for each algorithm
    // T_u of Dijkstra's algorithm is 0, T_q of Dijkstra's algorithm is dijk_query_times[i]
    // T_u of CGraph's algorithm is

    cout << "---------------------------------------------------------" << endl;
    cout << "uber's model:" << endl;

    for (int T_r = 100; T_r <= 10000; T_r *= 2) {
        cout << "T_r is " << T_r << endl;
        cout << "Object arrival rate: " << obj_num << " ";

        for (double T = 1.0; T <= 16.0; T *= 2) {
            cout << endl << "Period T is: " << T << endl;
            vector<double> Dijkstra_Thr, CGraph_Thr;
            double microsec_per_sec = 1000000.0;
            for (int alg = 0; alg < dijk_times.size(); alg++) {
                Dijkstra_Thr.push_back(
                        analytical_query_first(dijk_times[alg].query / microsec_per_sec, 0.0, T_r / microsec_per_sec,
                                               obj_num, T,
                                               dijk_times[alg].query_variance / microsec_per_sec / microsec_per_sec));
                CGraph_Thr.push_back(
                        analytical_query_first(cgraph_times[alg].query / microsec_per_sec,
                                               cgraph_times[alg].update / microsec_per_sec, T_r / microsec_per_sec,
                                               obj_num,
                                               T,
                                               cgraph_times[alg].query_variance / microsec_per_sec / microsec_per_sec));

            }
            double max_thr = -1;
            for (double thr : Dijkstra_Thr) {
                if (max_thr < thr) {
                    max_thr = thr;
                }
            }
            cout << "max throughput of Dijkstra's algorithm: " << max_thr << endl;

            double max_thr_ours = -1;
            int comb_id = -1;
            int ind = 0;
            for (double thr : CGraph_Thr) {
                if (max_thr_ours < thr) {
                    max_thr_ours = thr;
                    comb_id = ind;
                }
                ind++;
            }
            cout << "max throughput of CGraph's algorithm: " << max_thr_ours << endl;
            cout << "combination preferred is combination " << comb_id << endl;
            cout << endl;
            fprintf(file, "%f \t %f \t %d \t %f \t %f \t %d\n", obj_num / test_n, T, T_r, max_thr_ours, max_thr,
                    comb_id);
            accumulated_comb_ids.push_back(comb_id);
        }
    }
}

/**
 * Given the T_r and obj_ratio, output the simulation list of queries and updates
 * @param T_r
 * @param obj_ratio
 * @return
 */

// Using precomputed query times and update times
simulate_times test_uber_one_instance(int T_r, double obj_ratio, int k, double T) {
    simulate_times costs;

    cout << "---------------------------------------------------------" << endl;
    cout << "uber's model:" << endl;


    cout << "T_r is " << T_r << endl;
    cout << "Object arrival rate: " << test_n * obj_ratio << " ";
    string method_name = "ours";
    string train_input_file_name = make_out_file_name(method_name, k, obj_ratio, "training");
    string simulate_input_file_name = make_out_file_name(method_name, k, obj_ratio, "simulate");
    FILE *train_input_file = fopen((input_parameters.output_data_dir + train_input_file_name).c_str(), "r");
    if (train_input_file == NULL) {
        cout << train_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    FILE *simulate_input_file = fopen((input_parameters.output_data_dir + simulate_input_file_name).c_str(), "r");
    cout << "simulate input file: " << simulate_input_file_name << endl;
    if (simulate_input_file == NULL) {
        cout << simulate_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    fscanf(train_input_file, "%*[^\n]\n", NULL);
    int k_val;
    double obj_ratio_val;
    int comb_id_val;
    double query_val;
    double query_var_val;
    double update_val;
    double update_var_val;
    int max_comb_id = -1;
    double max_thr = -1.0;
    int cnt = 0;
    while (fscanf(train_input_file, "%d %lf %d %lf %lf %lf %lf", &k_val, &obj_ratio_val, &comb_id_val, &query_val,
                  &query_var_val,
                  &update_val, &update_var_val) == 7) {
        double thr = analytical_query_first(query_val * 1.0 / MICROSEC_PER_SEC, update_val * 1.0 / MICROSEC_PER_SEC,
                                            T_r * 1.0 / MICROSEC_PER_SEC,
                                            test_n * obj_ratio_val, T,
                                            query_var_val / MICROSEC_PER_SEC / MICROSEC_PER_SEC);
//        cout<<"thr: "<< thr<<endl;
        if (thr > max_thr) {
            max_thr = thr;
            max_comb_id = comb_id_val;
        }
//        cout<<++cnt<<endl;
        cout << k_val << " " << obj_ratio_val << " " << comb_id_val << " " << query_val << " " << query_var_val << " "
             << update_val << " " << update_var_val << endl;
    }
    cout << "max_thr: " << max_thr << endl;
    cout << "max id: " << max_comb_id << endl;
    stringstream ss;
    ss << max_comb_id;
    string started_line_str = "Combination " + ss.str() + "\n";
    char str[100];
    while (fgets(str, 100, simulate_input_file) != NULL) {
        string s(str);

        if (s.compare(started_line_str) != 0) continue;
        fgets(str, 100, simulate_input_file);
        cout << str << endl;
        while (fgets(str, 100, simulate_input_file) != NULL) {
//            string s(str);
//            if (s.compare("Update costs") != 0)
            if (str[0] <= '9' && str[0] >= '0') {
                int update_ind, cost;
                sscanf(str, "%d %d", &update_ind, &cost);
                while (costs.query_times.size() < update_ind + 1) {
                    vector<int> *new_vec = new vector<int>();
                    costs.query_times.push_back(*new_vec);
                }
                costs.query_times[update_ind].push_back(cost);

            } else {
//                cout<<"here2"<<str<<endl;
                while (fgets(str, 100, simulate_input_file) != NULL) {
//                    if (strncmp(str, "Combination", 10) != 0) {
//                    string s(str);
                    if (str[0] <= '9' && str[0] >= '0') {
                        costs.update_times.push_back(atoi(str));

                    } else {
                        break;
                    }

                }
                break;
            }

        }
        break;
    }
    fclose(train_input_file);
    fclose(simulate_input_file);
    return costs;

}

simulate_times test_query_first_by_scob_configuration(int durationTime, double obj_ratio, int k, double T,
                                                      int best_scob_configuration_id) {
    // Prepare
    prepare_environment();
    mem_struct mems;
    allocate_mem(mems, test_n + test_n + 1);
    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;
    int *object_nodes_hash = new int[test_n + 1];
    double obj_num = test_n * obj_ratio;

    // control test data
    obj_num = make_obj_data(obj_num, object_nodes_hash, objects_vector, 4);
    cout << (int) (obj_num) << " objects" << endl;
    vector<ScobConfiguration> configurations = ScobConfiguration::generate_all_configurations(level_num);

    // there are a number of (at most $index) algorithms, let's test them one by one

    tradeoff_level = configurations[best_scob_configuration_id].tradeoff_level;
    cut_level = configurations[best_scob_configuration_id].cut_level;
    mode = configurations[best_scob_configuration_id].mode;
    simulate_times simT;
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;

    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;

    int num_of_periods=(int) (durationTime / T);
    cout<<"num_of_periods: "<<num_of_periods<<endl;
    for (int period_id = 0; period_id <= num_of_periods; period_id++) {
        if(ZIPF_QUERY){
            if(period_id>0) break;
        }
//        cout<<period_id<<endl;
        /*
         * calculate the initial time, it is also the time spent on updating all the objects
         * each period_id (a duration of T), we initialize the kDNN once
         */
        if (test_parameters.data == DATA_BJ_REAL) {// for UCAR data
            if (test_parameters.arrival_model == QUEUE_POLICY_QUERY_FIRST) {
                obj_num = get_ucar_object_vector(car_period_location_list, object_nodes_hash, objects_vector,
                                                 period_id, indexes);

            }
        }
        vector<int> query_times_within_one_period;
        long start_init = clock();
        for (int i = 0; i <= test_n; i++) {
            hier_local_knn_arr[i].clear();
        }
        int is_init = 1;
        update_kdnn_list(is_init, k, 0, test_parameters.version, q, dist, visited, objects_vector);
        long scob_init_time = clock() - start_init;
        simT.update_times.push_back((int) (scob_init_time));
//        cout << "scob init time: " << scob_init_time << endl;
//        cout << "scob init time (per object): " << scob_init_time / obj_num << endl;

        // we don't consider the response requirement here, so that the collected queries would be more than enough
        double remain_time_within_one_period = (T * MICROSEC_PER_SEC - scob_init_time) * 1.0;
        // for each setting, we test "num" times, num by default = 200
        while (remain_time_within_one_period > 0) {
            int query=get_valid_random_node();
         //   cout<<"query: "<<query<<endl;
            long start = clock();
            scob_KNN_query(k, hier_local_knn_arr, query, dist, visited, q);
            long queryTime = clock() - start;
            remain_time_within_one_period -= queryTime;
            query_times_within_one_period.push_back((int) (queryTime));

        }
        simT.query_times.push_back(query_times_within_one_period);
    }
    delete_mems(mems);
    delete[] object_nodes_hash;
    return simT;
}



simulate_times
test_fifo_by_combination(int durationTime, double obj_ratio, int k, double update_obj_ratio, int max_comb_id) {
    // Prepare
    srand(time(NULL));
    int k_max = 50;
    simulate_times simT;
//    max_comb_id = 0;
    mem_struct mems;
    allocate_mem(mems, test_n + test_n + 1);
    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;
    int *car_nodes = new int[test_n + 1];
    double obj_num = test_n * obj_ratio;
    int query_num_per_update = 50;
    // control test data
    obj_num = make_obj_data(obj_num, car_nodes, objects_vector);
    cout << "the number of objects: " << objects_vector.size() << endl;
    vector<int> non_object_nodes;
    for (int i = 1; i < test_n; i++) {
        if (!car_nodes[i] && valid_node_hash[i])
            non_object_nodes.push_back(i);
    }
    //shuffle
//    for (int i = 0; i < 50000000; i++) {
//        // shuffle non objects
//        int z1 = rand() % non_object_nodes.size();
//        int z2 = rand() % non_object_nodes.size();
//        int u = non_object_nodes[z1];
//        non_object_nodes[z1] = non_object_nodes[z2];
//        non_object_nodes[z2] = u;
//
//        // shuffle objects
//        z1 = rand() % objects_vector.size();
//        z2 = rand() % objects_vector.size();
//        u = objects_vector[z1];
//        objects_vector[z1] = objects_vector[z2];
//        objects_vector[z2] = u;
//
//
//    }

    shuffle(objects_vector);
    shuffle(non_object_nodes);
    cout << (int) (obj_num) << " objects" << endl;
    int para1[100] = {0};
    int para2[100] = {0};
    int para3[100] = {0};
    int index = 0;

    // Comment the following a few lines, because the gentle t-climbs can be slow
    tradeoff_level = 0;
    // NOTE!!! Should comment for BJ
    if(network_name.compare("BJ")!=0) {
        for (cut_level = 0; cut_level <= level_num; cut_level++) {
            para1[index] = tradeoff_level;
            para2[index] = cut_level;
            para3[index++] = 2;
        }
    }
    // NOTE!!! Should comment for BJ
    tradeoff_level = level_num;
    for (cut_level = level_num; cut_level >= 0; cut_level--) {
        para1[index] = tradeoff_level;
        para2[index] = cut_level;
        para3[index++] = 3;
    }
    cout << "obj num: " << obj_num << endl;

    // there are a number of (at most $index) algorithms, let's test them one by one

    tradeoff_level = para1[max_comb_id];
    cut_level = para2[max_comb_id];
    mode = para3[max_comb_id];
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;


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

    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;


    for (int i = 0; i <= test_n; i++) {
        hier_local_knn_arr[i].clear();
    }
    update_kdnn_list(1, k_max, 0, test_parameters.version, q, dist, visited, objects_vector);
    // pokemon-go
    double remainTime = durationTime * 1000000;

    while (remainTime > 0) {
        vector<int> queryTimesPerUpdate;
        for (int z = 0; z < query_num_per_update; z++) {
            int query;
            do {

                query = rand() % test_n;
//                cout<<"query: "<<query<<endl;
            } while (valid_node_hash[query] == 0);
            long start = clock();
            vector<KNode> scobRes=scob_KNN_query(k, hier_local_knn_arr, query, dist, visited, q);

            if(verifyResult) {
                vector<KNode> naiveRes = naiveKNN(k, objects_vector, query, dist, visited, q, car_nodes);
                if(scobRes.size()!=naiveRes.size()){
                    cout<<"size different"<<endl;
                    stop_here();
                }
                for(int i = 0;i<scobRes.size();i++){
                    if(scobRes[i].dis!=naiveRes[i].dis){
                        cout<<"distance different"<<endl;
                        for(int j = 0;j<scobRes.size();j++) {
                            cout << j << " " << scobRes[j].id<<" "<<scobRes[j].dis << " " <<
                                 naiveRes[j].id<<" "<<naiveRes[j].dis << endl;
                        }
                        stop_here();
                    }
                }
            }
            remainTime -= (clock() - start);
            queryTimesPerUpdate.push_back((int) (clock() - start));


//            cout<<"query time: "<<clock()-start<<endl;
        }
        simT.query_times.push_back(queryTimesPerUpdate);
        long updateStart = clock();
        int pnum = rand() % 100;
        if ((pnum < 50 && (objects_vector.size() > 0) ) || (non_object_nodes.size() == 0) ) {
//            cout<<"remove start"<<endl;
//            cout<<"size: "<<objects_vector.size()<<endl;
            int removeNode = objects_vector[0];
//            cout<<"remove "<<removeNode<<endl;
            long start1 = clock();
            remove_kdnn_object(removeNode, k_max, k, mems, rev_shortcuts_arr, durationTime);
            remainTime -= clock() - start1;
            simT.update_times.push_back((int) (clock() - start1));
            remainTime -= (clock() - start1);
            objects_vector.erase(objects_vector.begin());
            car_nodes[removeNode] = 0;
            non_object_nodes.push_back(removeNode);
//            cout<<"remove end"<<endl;
        } else {
//            cout<<"insert start"<<endl;
            int inserted_node_id = non_object_nodes[0];
//            cout<<"insert "<<inserted_node_id<<endl;
            long start2 = clock();
            update_kDNN_object(0, inserted_node_id, k_max, k, dist, visited, q);
            remainTime -= clock() - start2;
            simT.update_times.push_back((int) (clock() - start2));
            remainTime -= (clock() - start2);
            non_object_nodes.erase(non_object_nodes.begin());
            objects_vector.push_back(inserted_node_id);
            car_nodes[inserted_node_id] = 1;
//            cout<<"insert end"<<endl;
        }


    }


    delete_mems(mems);
    delete[] car_nodes;
    cout << "finished simulations!" << endl;
    return simT;
}


simulate_times test_dijk_query_first_one_instance_online(int durationTime, int T_r, double obj_ratio, int k, double T) {
    cout << "---------------------------------------------------------" << endl;
    cout << "uber's model for Dijkstra:" << endl;


    cout << "T_r is " << T_r << endl;
    cout << "Object arrival rate: " << test_n * obj_ratio << " ";

    int isPOIInsert = 0;
    srand(time(NULL));
    read_road_network();
    mem_struct mems;
    allocate_mem(mems, test_n + test_n + 1);
    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;
    int *car_nodes = new int[test_n + 1];
    double obj_num = test_n * obj_ratio;

    // control test data
    if (test_parameters.data == 1) {// test shenzhou data
        obj_num = make_real_nodes("beijing_car_positions_node.txt", car_nodes, objects_vector);
    } else if (test_parameters.data == 2) {// test NW POI
        obj_num = make_real_nodes("NW_POI.txt", car_nodes, objects_vector);
        if (test_parameters.arrival_model == 0) {// poke model
            isPOIInsert = 1;
        }
    } else {
        obj_num = make_random_nodes((int) (obj_num), car_nodes, objects_vector);

    }
    cout << (int) (obj_num) << " objects" << endl;

    cout << "obj num: " << obj_num << endl;

    cout << "the number of objects: " << objects_vector.size() << endl;
    vector<int> non_object_nodes;
    if (isPOIInsert) {
        unordered_map<int, int> hash_nodes;
        for(int u: objects_vector) {
            int ind = rand()%level_scs_arr[0][u].size();

            int v = level_scs_arr[0][u][ind].id;
            if(!hash_nodes[v] && car_nodes[v]==0) {
                non_object_nodes.push_back(v);
                hash_nodes[v]=1;
            }

        }
    } else {
        for (int i = 1; i < test_n; i++) {
            if (!car_nodes[i] && valid_node_hash[i])
                non_object_nodes.push_back(i);
        }

    }




    //shuffle
//    for (int i = 0; i < 1000000; i++) {
//        // shuffle non objects
//        int z1 = rand() % non_object_nodes.size();
//        int z2 = rand() % non_object_nodes.size();
//        int u = non_object_nodes[z1];
//        non_object_nodes[z1] = non_object_nodes[z2];
//        non_object_nodes[z2] = u;
//
//        // shuffle objects
//        z1 = rand() % objects_vector.size();
//        z2 = rand() % objects_vector.size();
//        u = objects_vector[z1];
//        objects_vector[z1] = objects_vector[z2];
//        objects_vector[z2] = u;
//
//
//    }
    shuffle(non_object_nodes);
    shuffle(objects_vector);

    simulate_times simT;
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;


    for (int period = 0; period <= durationTime / T; period++) {
        /*
         * cal init time for C-Graph, it is also the time spent on updating all the objects
         */
        if(ZIPF_QUERY && period>0)
            break;
        vector<int> periodQueryTimes;
        long cgraph_init_time = 0;
        simT.update_times.push_back((int) (0));
//        cout << "cgraph init time: " << cgraph_init_time << endl;
//        cout << "cgraph init time (per object): " << cgraph_init_time / obj_num << endl;

        double remainTime = (T * 1000000 - cgraph_init_time) * 1.0;
        vector<int> valid_nodes;
        for(int i=0;i<test_n;i++){
            if(valid_node_hash[i]){
                valid_nodes.push_back(i);

            }
        }
//        cout<<"valid node size:"<<valid_nodes.size()<<endl;
//        char c;
//        cin>>c;
        // for each setting, we test "num" times, num by default = 200
        int flg=1;
        while (remainTime > 0) {
            int query;
            do {

                    query = valid_nodes[rand() % valid_nodes.size()];


            } while (valid_node_hash[query] == 0);
//            cout<<"query: "<<query<<endl;
            long start = clock();
            naiveKNN(k, objects_vector, query, dist, visited, q, car_nodes);
            long queryTime = clock() - start;

//            cout<<"query time: "<<queryTime<<endl;
            remainTime -= queryTime;
            periodQueryTimes.push_back((int) (queryTime));

            int pnum = rand() % 100;
            if ( (pnum < 50 && (objects_vector.size() > 0)) || (non_object_nodes.size() == 0)) {
//            if(flg){
                int removeNode = objects_vector[0];
                long start1 = clock();
                objects_vector.erase(objects_vector.begin());
//                cout<<"erase time1: "<<clock()-start1<<endl;
                car_nodes[removeNode] = 0;
                non_object_nodes.push_back(removeNode);
                flg=0;
            } else {
                int inserted_node_id = non_object_nodes[0];
                long start2 = clock();
                non_object_nodes.erase(non_object_nodes.begin());
//                cout<<"erase time2: "<<clock()-start2<<endl;
                objects_vector.push_back(inserted_node_id);
                car_nodes[inserted_node_id] = 1;
                flg=1;
            }
//            cout<<"obj num: "<<objects_vector.size()<<endl;

        }
        simT.query_times.push_back(periodQueryTimes);


        // estimate the updates of object nodes; since the cost is tiny, we need to test multiple times and
//        compute the average.
//        long updateStart = clock();
//        int tmp;
//        for (int updateAccess = 0; updateAccess < 10000; updateAccess++) {
//            tmp = car_nodes[updateAccess];
//        }
//        cout << tmp << endl;
//        long updateCost = (clock() - updateStart) / 10000.0;
//        cout << "update cost of Dijkstra: " << updateCost << endl;
//        simT.update_times.push_back(updateCost);
//        simT.update_times.push_back(0);

    }
    delete_mems(mems);
    delete[] car_nodes;
    return simT;


}



simulate_times test_dijk_fifo_one_instance_online(int durationTime, int T_r, double obj_ratio, int k, double T) {
    cout << "---------------------------------------------------------" << endl;
    cout << "poke's model for Dijkstra:" << endl;


    cout << "T_r is " << T_r << endl;
    cout << "Object arrival rate: " << test_n * obj_ratio << " ";

    int isPOIInsert = 0;
    srand(time(NULL));
    read_road_network();
    mem_struct mems;
    allocate_mem(mems, test_n + test_n + 1);
    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;
    int *car_nodes = new int[test_n + 1];
    double obj_num = test_n * obj_ratio;
    int query_num_per_update = 50;

    // control test data
    if (test_parameters.data == 1) {// test shenzhou data
        obj_num = make_real_nodes("beijing_car_positions_node.txt", car_nodes, objects_vector);
    } else if (test_parameters.data == 2) {// test NW POI
        obj_num = make_real_nodes("NW_POI.txt", car_nodes, objects_vector);
        if (test_parameters.arrival_model == 0) {// poke model
            isPOIInsert = 1;
        }
    } else {
        obj_num = make_random_nodes((int) (obj_num), car_nodes, objects_vector);

    }
    cout << (int) (obj_num) << " objects" << endl;

    cout << "obj num: " << obj_num << endl;

    cout << "the number of objects: " << objects_vector.size() << endl;
    vector<int> non_object_nodes;
    if (isPOIInsert) {
        unordered_map<int, int> hash_nodes;
        for(int u: objects_vector) {
            int ind = rand()%level_scs_arr[0][u].size();

            int v = level_scs_arr[0][u][ind].id;
            if(!hash_nodes[v] && car_nodes[v]==0) {
                non_object_nodes.push_back(v);
                hash_nodes[v]=1;
            }

        }
    } else {
        for (int i = 1; i < test_n; i++) {
            if (!car_nodes[i] && valid_node_hash[i])
                non_object_nodes.push_back(i);
        }

    }




    //shuffle
//    for (int i = 0; i < 1000000; i++) {
//        // shuffle non objects
//        int z1 = rand() % non_object_nodes.size();
//        int z2 = rand() % non_object_nodes.size();
//        int u = non_object_nodes[z1];
//        non_object_nodes[z1] = non_object_nodes[z2];
//        non_object_nodes[z2] = u;
//
//        // shuffle objects
//        z1 = rand() % objects_vector.size();
//        z2 = rand() % objects_vector.size();
//        u = objects_vector[z1];
//        objects_vector[z1] = objects_vector[z2];
//        objects_vector[z2] = u;
//
//
//    }
    shuffle(non_object_nodes);
    shuffle(objects_vector);

    simulate_times simT;
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;



    long cgraph_init_time = 0;
    simT.update_times.push_back((int) (0));
//        cout << "cgraph init time: " << cgraph_init_time << endl;
//        cout << "cgraph init time (per object): " << cgraph_init_time / obj_num << endl;

    double remainTime = (durationTime * 1000000) * 1.0;
    vector<int> valid_nodes;
    for(int i=0;i<test_n;i++){
        if(valid_node_hash[i]){
            valid_nodes.push_back(i);

        }
    }
//        cout<<"valid node size:"<<valid_nodes.size()<<endl;
//        char c;
//        cin>>c;
    // for each setting, we test "num" times, num by default = 200
    int flg=1;
    while (remainTime > 0) {
        vector<int> periodQueryTimes;
        for(int i = 0; i<query_num_per_update;i++) {
            int query;
            do {

                query = valid_nodes[rand() % valid_nodes.size()];


            } while (valid_node_hash[query] == 0);

//            cout<<"query: "<<query<<endl;
            long start = clock();
            naiveKNN(k, objects_vector, query, dist, visited, q, car_nodes);
            long queryTime = clock() - start;

//            cout<<"query time: "<<queryTime<<endl;
            remainTime -= queryTime;
            periodQueryTimes.push_back((int) (queryTime));
        }
        int pnum = rand() % 100;
        if ( (pnum < 50 && (objects_vector.size() > 0)) || (non_object_nodes.size() == 0)) {
//            if(flg){
            int removeNode = objects_vector[0];
            long start1 = clock();
            objects_vector.erase(objects_vector.begin());
//                cout<<"erase time1: "<<clock()-start1<<endl;
            car_nodes[removeNode] = 0;
            non_object_nodes.push_back(removeNode);
            flg=0;
        } else {
            int inserted_node_id = non_object_nodes[0];
            long start2 = clock();
            non_object_nodes.erase(non_object_nodes.begin());
//                cout<<"erase time2: "<<clock()-start2<<endl;
            objects_vector.push_back(inserted_node_id);
            car_nodes[inserted_node_id] = 1;
            flg=1;
        }
//            cout<<"obj num: "<<objects_vector.size()<<endl;
        simT.query_times.push_back(periodQueryTimes);
        simT.update_times.push_back(0);

    }
    delete_mems(mems);
    delete[] car_nodes;
    return simT;


}

int get_best_scob_configuration_by_training_data(FILE* train_input_file, int T_r, double T){
    fscanf(train_input_file, "%*[^\n]\n", NULL);
    int k_val;
    double obj_ratio_val;
    int comb_id_val;
    double query_val;
    double query_var_val;
    double update_val;
    double update_var_val;
    int best_scob_configuration_id = -1;
    double max_thr = -1.0;
    int cnt = 0;
    while (fscanf(train_input_file, "%d %lf %d %lf %lf %lf %lf", &k_val, &obj_ratio_val, &comb_id_val, &query_val,
                  &query_var_val,
                  &update_val, &update_var_val) == 7) {
        double thr = analytical_query_first(query_val * 1.0 / MICROSEC_PER_SEC, update_val * 1.0 / MICROSEC_PER_SEC,
                                            T_r * 1.0 / MICROSEC_PER_SEC,
                                            test_n * obj_ratio_val, T,
                                            query_var_val / MICROSEC_PER_SEC / MICROSEC_PER_SEC);
//        cout<<"thr: "<< thr<<endl;
        if (thr > max_thr) {
            max_thr = thr;
            best_scob_configuration_id = comb_id_val;
        }
//        cout<<++cnt<<endl;
        cout << k_val << " " << obj_ratio_val << " " << comb_id_val << " " << query_val << " " << query_var_val << " "
             << update_val << " " << update_var_val << endl;
    }
    cout << "max_thr: " << max_thr << endl;
    cout << "max id: " << best_scob_configuration_id << endl;
    return best_scob_configuration_id;
}
simulate_times test_toain_query_first_one_instance_online(int durationTime, int response_time_requirement_seconds,
                                                          double obj_ratio, int k, double period) {
    simulate_times costs;
    cout << "---------------------------------------------------------" << endl;
    cout << "Query first arrival model:" << endl;


    cout << "response time requirement is " << response_time_requirement_seconds << endl;
    cout << "Object arrival rate: " << test_n * obj_ratio << " ";
    string method_name = "ours";
    string train_input_file_name = make_out_file_name(method_name, k, obj_ratio, "training");
    FILE *train_input_file = fopen((input_parameters.output_data_dir + train_input_file_name).c_str(), "r");
    test_valid_file(train_input_file, (input_parameters.output_data_dir + train_input_file_name), "test_toain_query_first_one_instance_online");

    int best_scob_configuration_id = get_best_scob_configuration_by_training_data(train_input_file, response_time_requirement_seconds, period);

    // for coding simplicity, we separate the queueing simulation and online queries;
    // That is, we first collect all the necessary (actually,
    // more than enough) query
    // costs and update costs; and using those costs as input to the queuing simulator
    costs = test_query_first_by_scob_configuration(durationTime, obj_ratio, k, period, best_scob_configuration_id);

    fclose(train_input_file);

    return costs;

}




/**
 *
 * @param T_r
 * @param obj_ratio the ratio between the number of objects and the number of nodes
 * @param k
 * @param update_obj_ratio the ratio between the number of updated objects (per second) and the number of objects
 * @return
 */
simulate_times testing_poke_one_instance(int T_r, double obj_ratio, int k, double update_obj_ratio) {
    simulate_times costs;

    cout << "---------------------------------------------------------" << endl;
    cout << "poke's model:" << endl;


    cout << "T_r is " << T_r << endl;
    cout << "Total objects: " << test_n * obj_ratio << endl;
    cout << "Update arrival rate: " << test_n * obj_ratio / update_obj_ratio << " ";
    string method_name = "ours";
    string train_input_file_name = make_out_file_name(method_name, k, obj_ratio, "training");
    string simulate_input_file_name = make_out_file_name(method_name, k, obj_ratio, "simulate");
    FILE *train_input_file = fopen((input_parameters.output_data_dir + train_input_file_name).c_str(), "r");
    string test_output_file_name = make_out_file_name(method_name, k, obj_ratio, "testing");
    FILE *test_output_file = fopen((input_parameters.output_data_dir + test_output_file_name).c_str(), "w");
    if (test_output_file == NULL) {
        cout << test_output_file_name << " cannot open!" << endl;
        stop_here();
    }
    if (train_input_file == NULL) {
        cout << train_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    FILE *simulate_input_file = fopen((input_parameters.output_data_dir + simulate_input_file_name).c_str(), "r");
    if (simulate_input_file == NULL) {
        cout << simulate_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    cout << "open training file: " << train_input_file_name << endl;
    cout << "open simulating file: " << simulate_input_file_name << endl;
    fscanf(train_input_file, "%*[^\n]\n", NULL);
    int k_val;
    double obj_ratio_val;
    int comb_id_val;
    double query_val;
    double query_var_val;
    double update_val;
    double update_var_val;
    int max_comb_id = -1;
    double max_thr = -1.0;
    while (fscanf(train_input_file, "%d %lf %d %lf %lf %lf %lf", &k_val, &obj_ratio_val, &comb_id_val, &query_val,
                  &query_var_val,
                  &update_val, &update_var_val) == 7) {
        double thr = analytical_fifo(query_val * 1.0 / MICROSEC_PER_SEC, update_val * 1.0 / MICROSEC_PER_SEC,
                                     T_r * 1.0 / MICROSEC_PER_SEC,
                                     test_n * obj_ratio / update_obj_ratio,
                                     query_var_val / MICROSEC_PER_SEC / MICROSEC_PER_SEC,
                                     update_var_val / MICROSEC_PER_SEC / MICROSEC_PER_SEC);
//        cout<< thr<<endl;
        if (thr > max_thr) {
            max_thr = thr;
            max_comb_id = comb_id_val;
        }
    }

    cout << "max_comb_id: " << max_comb_id << " max_thr: " << max_thr << endl;
    stringstream ss;
    ss << max_comb_id;

    string started_line_str = "Combination " + ss.str() + "\n";
    char str[100];
    while (fgets(str, 100, simulate_input_file) != NULL) {
        string s(str);

        if (s.compare(started_line_str) != 0) continue;
        cout << str << endl;
        fgets(str, 100, simulate_input_file);

        while (fgets(str, 100, simulate_input_file) != NULL) {
//            string s(str);
//            if (s.compare("Update costs") != 0)
            if (str[0] <= '9' && str[0] >= '0') {
                int update_ind, cost;
                sscanf(str, "%d %d", &update_ind, &cost);
                while (costs.query_times.size() < update_ind + 1) {
                    vector<int> *new_vec = new vector<int>();
                    costs.query_times.push_back(*new_vec);
                }
                costs.query_times[update_ind].push_back(cost);

            } else {
                cout << "should be Update costs: " << str << endl;
                while (fgets(str, 100, simulate_input_file) != NULL) {
//                    if (strncmp(str, "Combination", 10) != 0) {
//                    string s(str);
                    if (str[0] <= '9' && str[0] >= '0') {
                        costs.update_times.push_back(atoi(str));

                    } else {
                        cout << "should be the next combination: " << str << endl;
                        break;
                    }

                }
                break;
            }

        }
        break;
    }
    int simulate_result = 0;
//    simulate_result=simulate(costs);
    fprintf(test_output_file, "%d %d\n", (int) (max_thr), simulate_result);
    fclose(train_input_file);
    fclose(simulate_input_file);
    fclose(test_output_file);
    return costs;

}

vector<simulate_times> testing_poke_one_instance_online(int T_r, double obj_ratio, int k, double update_obj_ratio) {
    vector<simulate_times> costs;
    int durationTime = 208;
    cout << "---------------------------------------------------------" << endl;
    cout << "poke's model:" << endl;

    cout << "T_r is " << T_r << endl;
    cout << "Total objects: " << test_n * obj_ratio << endl;
    cout << "Update arrival rate: " << test_n * obj_ratio / update_obj_ratio << " ";
    string method_name = "ours";
    string train_input_file_name = make_out_file_name(method_name, k, obj_ratio, "training");
    FILE *train_input_file = fopen((input_parameters.output_data_dir + train_input_file_name).c_str(), "r");
    if (train_input_file == NULL) {
        cout << train_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    cout << "open training file: " << train_input_file_name << endl;
    fscanf(train_input_file, "%*[^\n]\n", NULL);
    int k_val;
    double obj_ratio_val;
    int comb_id_val;
    double query_val;
    double query_var_val;
    double update_val;
    double update_var_val;
    int max_comb_id = -1;
    double max_thr = -1.0;
    int id_low=100, id_high=0;
    while (fscanf(train_input_file, "%d %lf %d %lf %lf %lf %lf", &k_val, &obj_ratio_val, &comb_id_val, &query_val,
                  &query_var_val,
                  &update_val, &update_var_val) == 7) {
        double thr = analytical_fifo(query_val * 1.0 / MICROSEC_PER_SEC, update_val * 1.0 / MICROSEC_PER_SEC,
                                     T_r * 1.0 / MICROSEC_PER_SEC,
                                     test_n * obj_ratio / update_obj_ratio,
                                     query_var_val / MICROSEC_PER_SEC / MICROSEC_PER_SEC,
                                     update_var_val / MICROSEC_PER_SEC / MICROSEC_PER_SEC);
        cout << thr << endl;
        if (thr > max_thr) {
            max_thr = thr;
            max_comb_id = comb_id_val;
        }
        if(id_low>comb_id_val) id_low = comb_id_val;
        if(id_high<comb_id_val) id_high = comb_id_val;
    }
//    max_comb_id=3;
    cout << "max_comb_id: " << max_comb_id << " max_thr: " << max_thr << endl;

    fclose(train_input_file);
//    max_comb_id = 7;
    int neighborhood_search_depth = 0;
    if(NEIGHBOR_SEARCH){
        neighborhood_search_depth=2;
    }

    cout<<"combinations on board: "<<endl;
    for(int comb_id = max(max_comb_id-neighborhood_search_depth, id_low);
        comb_id<=min(max_comb_id+neighborhood_search_depth, id_high); comb_id++) {
        cout<<"on boarding id "<<comb_id <<"... ";
        simulate_times temp_costs = test_fifo_by_combination(durationTime, obj_ratio, k, update_obj_ratio, comb_id);
        costs.push_back(temp_costs);
    }
    cout<<endl;

    return costs;

}



void test_toain_query_first() {
    int simulation_time = 208;
    test_parameters.arrival_model = QUEUE_POLICY_QUERY_FIRST;
    string method_name = "ours";
    int vary_k = 0;
    int vary_o = 0;
    int vary_T = 0;
    int vary_Tr = 1;
    string result_out_file_name;
    if (vary_k == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyk_result");
    if (vary_o == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyo_result");
    if (vary_T == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyT_result");
    if (vary_Tr == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyTr_result");
    if (vary_k == 0 && vary_o == 0 && vary_T == 0 && vary_Tr == 0)
        result_out_file_name = "case_study";
    FILE *result_out_file = fopen((input_parameters.output_data_dir + result_out_file_name).c_str(), "w");
    test_valid_file(result_out_file, (input_parameters.output_data_dir + result_out_file_name), "test_toain_query_first");
    fprintf(result_out_file, "T_r \t T \t k \t \t Obj_ratio \t Thr\n");

    vector<double> obj_ratio_vals = {0.001, 0.01, 0.1};
//    vector<double> obj_ratio_vals_def = {3000.0/test_n};
//    vector<double> obj_ratio_vals_def = {15000.0 / test_n};// case 2
    vector<double> obj_ratio_vals_def = {0.001};
//    vector<double> obj_ratio_vals_def = {5000.0 / test_n};
    vector<int> k_vals = {1, 10, 20, 30, 40};
    vector<int> k_vals_def = {20};
//    vector<int> k_vals_def = {1};
    vector<int> Tr_vals = {100, 400, 1600, 6400};
    vector<int> Tr_vals_def = {800};
    vector<double> T_vals = {1.0, 4.0, 16.0};
    vector<double> T_vals_def = {4.0};

    vector<double> selected_objs = obj_ratio_vals_def, selected_T_vals = T_vals_def;
    vector<int> selected_k_vals = k_vals_def, selected_Tr_vals = Tr_vals_def;
    if (vary_o) selected_objs = obj_ratio_vals;
    if (vary_k) selected_k_vals = k_vals;
    if (vary_T) selected_T_vals = T_vals;
    if (vary_Tr) selected_Tr_vals = Tr_vals;

    int tests = 1;
    vector<double> ths;
    for (double obj_ratio:selected_objs) {
        for (int k:selected_k_vals) {
            for (int T_r : selected_Tr_vals) {
                for (double T: selected_T_vals) {
                    for(int test_i=0;test_i<tests;test_i++) {
                        simulate_times selected_result = test_toain_query_first_one_instance_online(simulation_time,
                                                                                                    T_r,
                                                                                                    obj_ratio,
                                                                                                    k, T);
//                    for(int max_comb_id =0;max_comb_id<=13;max_comb_id+=3) {
//                        simulate_times selected_result = test_query_first_by_scob_configuration(simulation_time, obj_ratio, k, T,
//                                                                                max_comb_id);

                        cout << "finish online test..." << endl;
                        vector<double> simulated_query_costs;
                        vector<double> simulated_update_costs;
                        int th = test_for_queryfirst_divided_by_updates(simulated_query_costs, simulated_update_costs,
                                                                        simulation_time, T_r, T * MICROSEC_PER_SEC,
                                                                        selected_result.query_times,
                                                                        selected_result.update_times);
                        ths.push_back(th*1.0);
                        if (test_parameters.data == DATA_BJ_REAL) obj_ratio = BJ_OBJECTS * 1.0 / test_n;
                        fprintf(result_out_file, "%d \t %lf \t %d \t %lf \t %d\n", T_r, T, k, obj_ratio, th);
                        cout << "query mean: " << get_mean_double(simulated_query_costs) << endl;
                        cout << "update mean: " << get_mean_double(simulated_update_costs) << endl;
                        cout << "query variance: " << get_var_double(simulated_query_costs) << endl;
                        cout << "update variance: " << get_var_double(simulated_update_costs) << endl;

                        cout << T_r << " " << T << " " << k << " " << obj_ratio << " " << th << endl;
//                    }
                    }
                    cout<<"standard variance: "<<sqrt(get_var_double(ths))<<endl;
                }
            }
        }
    }
    fclose(result_out_file);
}



void test_dijk_query_first() {
    int simulation_time = 208;
    prepare_environment();
    test_parameters.arrival_model = QUEUE_POLICY_QUERY_FIRST;
    string method_name = "dijk";

    int vary_k = 0;
    int vary_o = 0;
    int vary_T = 0;
    int vary_Tr = 1;
    string result_out_file_name;
    if (vary_k == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyk_result");
    if (vary_o == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyo_result");
    if (vary_T == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyT_result");
    if (vary_Tr == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyTr_result");
    if (vary_k == 0 && vary_o == 0 && vary_T == 0 && vary_Tr == 0)
        result_out_file_name = "case_study";
    FILE *result_out_file = fopen((input_parameters.output_data_dir + result_out_file_name).c_str(), "w");
    if (result_out_file == NULL) {
        cout << "cannot open file: " << result_out_file_name << " in test_toain_query_first()" << endl;
        stop_here();
    }
    fprintf(result_out_file, "T_r \t T \t k \t \t Obj_ratio \t Thr\n");

    vector<double> obj_ratio_vals = {0.001, 0.01, 0.1};
    vector<double> obj_ratio_vals_def = {0.001};
//    vector<double> obj_ratio_vals_def = {15000.0 / test_n};
//    vector<double> obj_ratio_vals_def = {5000.0 / test_n};
    vector<int> k_vals = {1, 10, 20, 30, 40};
//    vector<int> k_vals_def = {1};
    vector<int> k_vals_def = {20};
//    vector<int> Tr_vals = {100, 200, 400, 800, 1600, 3200, 6400};
    vector<int> Tr_vals = {3200, 6400};
    vector<int> Tr_vals_def = {800};
    vector<double> T_vals = {1.0, 2.0, 4.0, 8.0, 16.0};
//    vector<double> T_vals_def = {4.0};
    vector<double> T_vals_def = {4.0};

    vector<double> selected_objs = obj_ratio_vals_def, selected_T_vals = T_vals_def;
    vector<int> selected_k_vals = k_vals_def, selected_Tr_vals = Tr_vals_def;
    if (vary_o) selected_objs = obj_ratio_vals;
    if (vary_k) selected_k_vals = k_vals;
    if (vary_T) selected_T_vals = T_vals;
    if (vary_Tr) selected_Tr_vals = Tr_vals;

    for (double obj_ratio:selected_objs) {
        for (int k:selected_k_vals) {
            for (int T_r : selected_Tr_vals) {
                for (double T: selected_T_vals) {

//                    simulate_times selected_result = test_uber_one_instance(T_r, obj_ratio, k, T);
                    simulate_times selected_result = test_dijk_query_first_one_instance_online(simulation_time, T_r, obj_ratio,
                                                                                   k, T);
                    vector<double> simulated_query_costs;
                    vector<double> simulated_update_costs;
                    int th = test_for_queryfirst_divided_by_updates(simulated_query_costs, simulated_update_costs,
                                                                    simulation_time, T_r, T * MICROSEC_PER_SEC,
                                                                    selected_result.query_times,
                                                                    selected_result.update_times);
                    if (test_parameters.data == 1) obj_ratio = BJ_OBJECTS * 1.0 / test_n;
                    fprintf(result_out_file, "%d \t %lf \t %d \t %lf \t %d\n", T_r, T, k, obj_ratio, th);

                    cout << "query mean: " << get_mean_double(simulated_query_costs) << endl <<
                         " update mean: " << get_mean_double(simulated_update_costs) << endl
                         << "query variance: " << get_var_double(simulated_query_costs)
                         << "update variance: " << get_var_double(simulated_update_costs) << endl;
                    cout << T_r << " " << T << " " << k << " " << obj_ratio << " " << th << endl;

                }
            }
        }
    }
    fclose(result_out_file);
}


void model_correctness_query_first() {
    test_parameters.arrival_model = 1;
    string method_name = "ours";

    for (double obj_ratio = 0.001; obj_ratio <= 0.001; obj_ratio *= 10) {
//        vector<int> k_vals = {1, 10, 20, 30, 40};
        vector<int> k_vals = {20};
        for (int k:k_vals) {
            for (int T_r = 800; T_r <= 800; T_r *= 2) {
                for (double T = 1.0; T <= 16.0; T *= 2) {
                    string result_out_file_name = make_out_file_name(method_name, k, obj_ratio, "result");
                    FILE *result_out_file = fopen((input_parameters.output_data_dir + result_out_file_name).c_str(), "w");
                    if (result_out_file == NULL) {
                        cout << "cannot open file: " << result_out_file_name << " in test_toain_query_first()" << endl;
                        stop_here();
                    }
                    simulate_times selected_result = test_uber_one_instance(T_r, obj_ratio, k, T);

                    vector<double> simulated_query_costs;
                    vector<double> simulated_update_costs;
                    int th = test_for_queryfirst_divided_by_updates(simulated_query_costs, simulated_update_costs,
                                                                    T, T_r, T * MICROSEC_PER_SEC,
                                                                    selected_result.query_times,
                                                                    selected_result.update_times);
                    cout << T_r << " " << T << " " << k << " ";
                    cout << th << endl;

                    double t_q = get_mean_double(simulated_query_costs);

                    double v_q = get_var_double(simulated_query_costs);

                    double t_u = simulated_update_costs[0];

                    cout << "update size: " << simulated_update_costs.size() << endl;

                    cout << t_q << " " << v_q << " " << t_u << endl;

                    if (test_parameters.data == 1) obj_ratio = 2770 * 1.0 / test_n;

                    double f1_thr = analytical_query_first(t_q * 1.0 / MICROSEC_PER_SEC,
                                                           t_u / (test_n * obj_ratio) / MICROSEC_PER_SEC,
                                                           T_r * 1.0 / MICROSEC_PER_SEC,
                                                           test_n * obj_ratio, T,
                                                           v_q / MICROSEC_PER_SEC / MICROSEC_PER_SEC);
                    cout << "f1_thr: " << f1_thr << " " << (f1_thr - th) / th << endl;


                    fclose(result_out_file);
                }
            }
        }
    }
}


void test_toain_fifo() {
    prepare_environment();
    test_parameters.arrival_model = QUEUE_POLICY_FIFO;
    string method_name = "ours";
    int vary_k = 0;
    int vary_o = 1;
    int vary_update_ratio = 0;
    int vary_Tr = 0;
    string result_out_file_name;
    if (vary_k == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyk_result");
    if (vary_o == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyo_result");
    if (vary_update_ratio == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyUr_result");
    if (vary_Tr == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyTr_result");
    if (vary_k == 0 && vary_o == 0 && vary_update_ratio == 0 && vary_Tr == 0) {
        result_out_file_name = "case_study";
    }
    FILE *result_out_file = fopen((input_parameters.output_data_dir + result_out_file_name).c_str(), "w");
    if (result_out_file == NULL) {
        cout << "cannot open file: " << result_out_file_name << " in test_toain_fifo()" << endl;
        stop_here();
    }
    fprintf(result_out_file, "T_r \t update_arrival \t k \t \t Obj_ratio \t Thr\n");

    if(vary_Tr || vary_update_ratio) NEIGHBOR_SEARCH=1;
    else NEIGHBOR_SEARCH=0;
    vector<double> obj_ratio_vals = {0.1, 0.01, 0.001};
    vector<double> obj_ratio_vals_def = {0.001};
//    vector<double> obj_ratio_vals_def = {5000.0/test_n};
    vector<int> k_vals = {1, 10, 20, 30, 40};
    vector<int> k_vals_def = {20};
//    vector<int> k_vals_def = {9};
    vector<int> Tr_vals = {6400, 3200, 1600, 800, 400, 200};
    vector<int> Tr_vals_def = {800};
    vector<double> Ur_vals = {1.0/16, 2.0/16, 4.0/16, 8.0/16, 1.0, 2.0, 4.0, 8.0, 16.0};
//    vector<double> Ur_vals = {8.0, 16.0};
    vector<double> Ur_vals_def = {4};
//    vector<double> Ur_vals_def = {5000.0/50000.0};

    vector<double> selected_objs = obj_ratio_vals_def, selected_Ur_vals = Ur_vals_def;
    vector<int> selected_k_vals = k_vals_def, selected_Tr_vals = Tr_vals_def;
    if (vary_o) selected_objs = obj_ratio_vals;
    if (vary_k) selected_k_vals = k_vals;
    if (vary_update_ratio) selected_Ur_vals = Ur_vals;
    if (vary_Tr) selected_Tr_vals = Tr_vals;
    int tests=1;
    for (double obj_ratio:selected_objs) {
        for (int k:selected_k_vals) {
            for (int T_r : selected_Tr_vals) {
                for (double update_obj_ratio: selected_Ur_vals) {
                    int total_th=0;
                    for(int i=0;i<tests;i++) {
                        if (test_parameters.data == 1) obj_ratio = 2770 * 1.0 / test_n;// BJ_real
                        if (test_parameters.data == 2) obj_ratio = 13000 * 1.0 / test_n;// NW_real

                        vector<simulate_times> selected_result_list = testing_poke_one_instance_online(T_r, obj_ratio, k,
                                                                                          update_obj_ratio);

                        int max_th = 0;

                        for(simulate_times selected_result : selected_result_list) {
                            cout << selected_result.query_times.size() << " " << selected_result.update_times.size()
                                 << endl;
                            vector<double> simulated_query_costs;
                            vector<double> simulated_update_costs;

                            if (selected_result.query_times.size() != 0 && selected_result.update_times.size() != 0) {
                                int th = test_for_fifo_divided_by_updates(simulated_query_costs, simulated_update_costs,
                                                                      208, T_r, test_n * obj_ratio / update_obj_ratio,
                                                                      selected_result.query_times,
                                                                      selected_result.update_times);
                                cout << " throughput: "<< th << endl << "query mean: " << get_mean_double(simulated_query_costs) << endl <<
                                     " update mean: " << get_mean_double(simulated_update_costs) << endl <<
                                     "query size: " << simulated_query_costs.size() << endl <<
                                     "update size:" << simulated_update_costs.size() << endl;

                                if(th > max_th) max_th = th;
//
//                            th = test_for_fifo_divided_by_updates(simulated_query_costs, simulated_update_costs,
//                                                                  128, T_r, test_n * obj_ratio / update_obj_ratio*10,
//                                                                  selected_result.query_times,
//                                                                  selected_result.update_times);
//                            cout<<th<<endl;
//                            cout << "query mean: " << get_mean_double(simulated_query_costs) << endl <<
//                                 " update mean: " << get_mean_double(simulated_update_costs) << endl;
//                            char c;
//                            cin>>c;
                            }

//                        cout << "query mean: " << get_mean_double(simulated_query_costs) << endl <<
//                             " update mean: " << get_mean_double(simulated_update_costs) << endl <<
//                             "query size: "<<simulated_query_costs.size()<<endl<<
//                             "update size:"<<simulated_update_costs.size()<<endl;
//                            fprintf(result_out_file, "%d \t %lf \t %d \t %lf \t %d\n", T_r, update_obj_ratio, k,
//                                    obj_ratio,
//                                    th);
//                        if(th==0) break;// no need to tests more
                            selected_result.query_times.clear();
                            selected_result.update_times.clear();
                            simulated_query_costs.clear();
                            simulated_update_costs.clear();
                        }
                        total_th += max_th;
                    }

                        cout << T_r << " " << update_obj_ratio << " " << k << " " << obj_ratio << " " << total_th/tests << endl;


                }
            }
        }
    }
    fclose(result_out_file);
}


void test_dijk_fifo() {
    cout<<"start test_dijk_fifo()"<<endl;
    int simulation_time = 208;
    prepare_environment();

    test_parameters.arrival_model = QUEUE_POLICY_FIFO;
    string method_name = "dijk";
    int vary_k = 0;
    int vary_o = 0;
    int vary_update_ratio = 0;
    int vary_Tr = 1;
    string result_out_file_name;
    if (vary_k == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyk_result");
    if (vary_o == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyo_result");
    if (vary_update_ratio == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyUr_result");
    if (vary_Tr == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyTr_result");
    if (vary_k == 0 && vary_o == 0 && vary_update_ratio == 0 && vary_Tr == 0)
        result_out_file_name = "case_study";
    FILE *result_out_file = fopen((input_parameters.output_data_dir + result_out_file_name).c_str(), "w");
    if (result_out_file == NULL) {
        cout << "cannot open file: " << result_out_file_name << " in test_toain_fifo()" << endl;
        stop_here();
    }
    fprintf(result_out_file, "T_r \t update_arrival \t k \t \t Obj_ratio \t Thr\n");

    vector<double> obj_ratio_vals = {0.1, 0.01, 0.001};
//    vector<double> obj_ratio_vals_def = {200000.0/test_n};
//    vector<double> obj_ratio_vals_def = {13000.0 / test_n};
//    vector<double> obj_ratio_vals_def = {5000.0 / test_n};
    vector<double> obj_ratio_vals_def = {0.001};
    vector<int> k_vals = {1, 10, 20, 30, 40};
//    vector<int> k_vals_def = {9};// case 1
    vector<int> k_vals_def = {20};
//    vector<int> Tr_vals = {100, 200, 400, 800, 1600, 3200, 6400};
    vector<int> Tr_vals = {3200, 6400};
    vector<int> Tr_vals_def = {800};
    vector<double> Ur_vals = {1.0, 2.0, 4.0, 8.0, 16.0};
    vector<double> Ur_vals_def = {4.0};
//    vector<double> Ur_vals_def = {5000.0/100000};

    vector<double> selected_objs = obj_ratio_vals_def, selected_Ur_vals = Ur_vals_def;
    vector<int> selected_k_vals = k_vals_def, selected_Tr_vals = Tr_vals_def;
    if (vary_o) selected_objs = obj_ratio_vals;
    if (vary_k) selected_k_vals = k_vals;
    if (vary_update_ratio) selected_Ur_vals = Ur_vals;
    if (vary_Tr) selected_Tr_vals = Tr_vals;

    int test_num = 1;
    for (double obj_ratio:selected_objs) {
        for (int k:selected_k_vals) {
            for (int T_r : selected_Tr_vals) {
                for (double update_obj_ratio: selected_Ur_vals) {
//                    cout<<"here"<<endl;
                    long total_th=0;
                    for(int test=0;test<test_num;test++) {
                        if (test_parameters.data == DATA_BJ_REAL) obj_ratio = BJ_OBJECTS * 1.0 / test_n;// BJ_real
                        if (test_parameters.data == DATA_NW_REAL) obj_ratio = NW_OBJECTS * 1.0 / test_n;// NW_real

                        simulate_times selected_result = test_dijk_fifo_one_instance_online(simulation_time, T_r, obj_ratio,
                                                                                       k, update_obj_ratio);

                        cout << "query size: " << selected_result.query_times.size() << " update size: "
                             << selected_result.update_times.size() << endl;
//                    for(int zz=0;zz<selected_result.update_times.size();zz++){
//                        cout<<"update: "<< selected_result.update_times[zz]<<endl;
//                        cout<<"query: "<<selected_result.query_times[zz][0]<<endl;
//                    }
                        vector<double> simulated_query_costs;
                        vector<double> simulated_update_costs;
                        int th = 0;
                        if (selected_result.query_times.size() != 0 && selected_result.update_times.size() != 0)
                            th = test_for_fifo_divided_by_updates(simulated_query_costs, simulated_update_costs,
                                                                  simulation_time, T_r,
                                                                  test_n * obj_ratio / update_obj_ratio,
                                                                  selected_result.query_times,
                                                                  selected_result.update_times);
                        total_th+=th;
                        cout << "query mean: " << get_mean_double(simulated_query_costs) << endl <<
                             " update mean: " << get_mean_double(simulated_update_costs) << endl
                             << "query variance: " << get_var_double(simulated_query_costs)
                             << "update variance: " << get_var_double(simulated_update_costs) << endl;
                        cout << T_r << " " << update_obj_ratio << " " << k << " " << obj_ratio << " " << th << endl;


                    }
                    cout << ALPHA<<" "<<T_r << " " << update_obj_ratio << " " << k << " " << obj_ratio << " " << total_th/test_num << endl;

//                    fprintf(result_out_file, "%d \t %lf \t %d \t %lf \t %d\n", T_r, update_obj_ratio, k, obj_ratio,
//                            th);



                }
            }
        }
    }
    fclose(result_out_file);
}





void model_correct_poke() {
    srand(time(NULL));
    test_parameters.arrival_model = 0;
    string method_name = "ours";

    for (double obj_ratio = 0.001; obj_ratio <= 0.001; obj_ratio *= 10) {
//        vector<int> k_vals = {1, 10, 20, 30, 40};
        vector<int> k_vals = {20};
        for (int k:k_vals) {
            for (int T_r = 800; T_r <= 800; T_r *= 2) {
                for (double Ur = 4.0; Ur <= 4.0; Ur *= 2) {
//                    string result_out_file_name = make_out_file_name(method_name, k, obj_ratio, "result");
//                    FILE* result_out_file = fopen((output_data_dir+result_out_file_name).c_str(), "w");
//                    if(result_out_file==NULL){
//                        cout<<"cannot open file: "<<result_out_file_name<<" in model_correct_poke()"<<endl;
//                        stop_here();
//                    }
                    if (test_parameters.data == 1) obj_ratio = 2770 * 1.0 / test_n;

                    simulate_times selected_result = testing_poke_one_instance(T_r, obj_ratio, k, Ur);
//                    simulate_times selected_result;
//                    for(int zz=0;zz<100;zz++){
//                        selected_result.update_times.push_back(1000);
//                        vector<int> test_vec;
//                        for(int hh=0;hh<100;hh++)
//                            test_vec.push_back(700);
//                        selected_result.query_times.push_back(test_vec);
//
//                    }

                    vector<double> simulated_query_costs;
                    vector<double> simulated_update_costs;
                    cout << "selected query size: " << selected_result.query_times.size() << endl;
                    cout << "selected update size: " << selected_result.update_times.size() << endl;
                    cout << "two means: " << get_mean_int(selected_result.query_times[0]) << " "
                         << get_mean_int(selected_result.update_times)
                         << endl;
                    cout << "object update ratio: " << test_n * obj_ratio / Ur << endl;
                    double th_all = 0.0;
                    for (int test_num = 0; test_num < 30; test_num++) {
                        th_all += test_for_fifo_divided_by_updates(simulated_query_costs, simulated_update_costs,
                                                                   1, T_r, test_n * obj_ratio / Ur,
                                                                   selected_result.query_times,
                                                                   selected_result.update_times);
                    }
                    cout << T_r << " " << Ur << " " << k << " ";
                    double th = th_all / 30.0;
                    cout << th << endl;

                    cout << "simulated_query_costs.size(): " << simulated_query_costs.size() << endl;
                    cout << "simulated_update_costs.size(): " << simulated_update_costs.size() << endl;
                    double t_q = get_mean_double(simulated_query_costs);

                    double v_q = get_var_double(simulated_query_costs);

                    double t_u = get_mean_double(simulated_update_costs);

                    double v_u = get_var_double(simulated_update_costs);

                    cout << "update size: " << simulated_update_costs.size() << endl;

                    cout << t_q << " " << v_q << " " << t_u << " " << v_u << endl;

                    double f1_thr = analytical_fifo(t_q * 1.0 / MICROSEC_PER_SEC, t_u / MICROSEC_PER_SEC,
                                                    T_r * 1.0 / MICROSEC_PER_SEC,
                                                    test_n * obj_ratio / Ur,
                                                    v_q / MICROSEC_PER_SEC / MICROSEC_PER_SEC,
                                                    v_u / MICROSEC_PER_SEC / MICROSEC_PER_SEC);

                    cout << "f1_thr: " << f1_thr << " " << 1 - (f1_thr - th) / th << endl;


//                    fclose(result_out_file);
                }
            }
        }
    }
}

//
//void test_toain_fifo() {
//    test_parameters.arrival_model=0;
//    string method_name = "ours";
//
//    for (double obj_ratio = 0.001; obj_ratio <= 0.001; obj_ratio *= 10) {
//        vector<int> k_vals = {1, 10, 20, 30, 40};
//        for (int k:k_vals) {
//            for (int T_r = 100; T_r <= 10000; T_r *= 2) {
//                for (double update_obj_ratio = 1.0; update_obj_ratio <= 16.0; update_obj_ratio *= 2) {
//                    string result_out_file_name = make_out_file_name(method_name, k, obj_ratio, "result");
//                    FILE* result_out_file = fopen((output_data_dir+result_out_file_name).c_str(), "w");
//                    if(result_out_file==NULL){
//                        cout<<"cannot open file: "<<result_out_file_name<<" in test_toain_query_first()"<<endl;
//                        stop_here();
//                    }
//                    simulate_times selected_result = testing_poke_one_instance(T_r, obj_ratio, k, update_obj_ratio);
////                    fprintf(result_out_file, "Query costs\n");
////                    for(int query_cost: selected_result.query_times)
////                        fprintf(result_out_file, "%d\n", query_cost);
////                    fprintf(result_out_file, "Update costs\n");
////                    for(int update_cost: selected_result.update_times)
////                        fprintf(result_out_file, "%d\n", update_cost);
//                    fclose(result_out_file);
//                }
//            }
//        }
//    }
//}

#endif //SOB_TESTING_H
