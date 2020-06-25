//
// Created by lostrong on 17/3/8.
//

#ifndef SOB_SIMULATION_H
#define SOB_SIMULATION_H

#include "../para/TestParameter.h"
#include "../para/GlobalVariables.h"
#include "../RoadKNNUtilities.h"
#include "../scob/ScobIndex.h"
#include "../util/DataFunctions.h"
#include "../util/TimeStat.h"
#include "../scob/ScobQuery.h"


simulate_times collect_simulate_poke_one_instance(int k_max, int k, int *&car_nodes,
                                                  mem_struct &mems, double obj_num, int num_of_updates,
                                                  int num_of_queries_per_update, int simulation_duration) {
    int  simulation_duration_micro_secs= simulation_duration * MICROSEC_PER_SEC;
    simulate_times cgraph_current_times;
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;
    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;
    long start_init = clock();
    for (int i = 0; i <= test_n; i++) {
        hier_local_knn_arr[i].clear();
    }

    update_kdnn_list(1, k_max, 0, test_parameters.version, q, dist, visited, objects_vector);

    long cgraph_init_time = clock() - start_init;
    cout << "cgraph init time: " << cgraph_init_time << endl;
    long duration = 0;


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

    for (int i = 0; i < num_of_updates; i++) {

        long duration_start=clock();
        int action = 0;
        int affected_node;
        long start2 = clock();

        int pnum = rand()%100;
        if (pnum<50) {
            remove_kdnn_object(objects_vector[0], k_max, k, mems, rev_shortcuts_arr, 10);
            action = -1;            //remove
            affected_node = objects_vector[0];

        } else {
            int inserted_node_id;
            do {
                inserted_node_id = get_rand_node();
            } while (car_nodes[inserted_node_id] || !valid_node_hash[inserted_node_id]);
            update_kDNN_object(0, inserted_node_id, k_max, k, dist, visited, q);
            action = 1;
            affected_node = inserted_node_id;

        }


        cgraph_current_times.update_times.push_back(clock() - start2);


        if (action == 1) {
            objects_vector.push_back(affected_node);
            car_nodes[affected_node] = 1;
        }
        if (action == -1) {
            objects_vector.erase(objects_vector.begin());
            car_nodes[affected_node] = 0;
        }



        /*
         * Method: C-Graph; Purpose: estimate the query time
         */
        vector<int> query_costs_per_update;
        for (int j = 0; j < num_of_queries_per_update; j++) {

            int query=get_valid_random_node();
            long start = clock();
            vector<KNode> r2;
            r2 = scob_KNN_query(k, hier_local_knn_arr, query, dist, visited, q);
            query_costs_per_update.push_back(clock() - start);

        }
        cgraph_current_times.query_times.push_back(query_costs_per_update);
        duration+=clock()-duration_start;
        if(duration> simulation_duration_micro_secs)
            break;

    }
    cgraph_current_times.k = k;
    cgraph_current_times.network_name = network_name;
    if (test_parameters.data == 0)
        cgraph_current_times.obj_name = "ran";
    else
        cgraph_current_times.obj_name = "real";
    cgraph_current_times.method_name = "ours";
    cgraph_current_times.obj_ratio = process(obj_num * 1.0 / test_n, 6);


//	cal_uber_stat_only(cgraph_times, obj_num, accumulated_comb_ids);
    return cgraph_current_times;
}


simulate_times collect_simulate_uber_one_instance(int *&car_nodes, int k,
                                                  mem_struct &mems, double obj_num, int num_of_periods,
                                                  int min_avg_queries_per_period, int duration_per_period) {
    int  duration_limit_per_period= duration_per_period * MICROSEC_PER_SEC;
    simulate_times cgraph_current_times;
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;
    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;
    /*
     * cal init time for C-Graph, it is also the time spent on updating all the objects
     */


    for (int pr = 0; pr < num_of_periods; pr++) {
        long duration = 0;

        make_obj_data(obj_num, car_nodes, objects_vector);
        long start_init = clock();

        for (int i = 0; i <= test_n; i++) {
            hier_local_knn_arr[i].clear();
        }
        update_kdnn_list(1, k, 0, test_parameters.version, q, dist, visited, objects_vector);
        long cgraph_init_time = clock() - start_init;
        cout << "cgraph init time: " << cgraph_init_time << endl;
        cout << "cgraph init time (per object): " << cgraph_init_time / obj_num << endl;
        cgraph_current_times.update_times.push_back(cgraph_init_time);// batch updates;
        duration+=cgraph_init_time;
        vector<int> query_times_per_update;
        int qr=0;
        while(true) {

            if(qr>min_avg_queries_per_period && duration > duration_limit_per_period)
                break;
                // for each setting, we test "num" times, num by default = 200
                int query = get_rand_node();

                /*
                 * Method: C-Graph; Purpose: estimate the query time
                 */
                long start = clock();
                vector<KNode> r2;
                r2 = scob_KNN_query(k, hier_local_knn_arr, query, dist, visited, q);

                query_times_per_update.push_back(clock() - start);
            duration+=clock()-start;
            qr++;

        }
        cgraph_current_times.query_times.push_back(query_times_per_update);
    }

    cgraph_current_times.k = k;
    cgraph_current_times.network_name = network_name;
    if (test_parameters.data == 0)
        cgraph_current_times.obj_name = "ran";
    else
        cgraph_current_times.obj_name = "real";
    cgraph_current_times.method_name = "ours";
    cgraph_current_times.obj_ratio = process(obj_num * 1.0 / test_n, 6);


//	cal_uber_stat_only(cgraph_times, obj_num, accumulated_comb_ids);
    return cgraph_current_times;
}


void simulate_toain_query_first() {
    // Prepare
    test_parameters.arrival_model = 1;
    srand(time(NULL));
    read_road_network();
    set_base_cell_num();
    int max_taxi_num = test_n;
    read_scob_index();
    mem_struct mems;
    allocate_mem(mems, test_n + 1 + max_taxi_num);

    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;

    int *car_nodes = new int[test_n + 1];
    vector<double> num_objs;

    if (test_parameters.data == 0)
        for (double nums = test_n * 0.001; nums <= test_n * 0.001; nums *= 10)
            num_objs.push_back(nums);
    else
        num_objs.push_back(0);

    int num_of_periods = 10;
    int min_avg_queries_per_period = 20;
    for (double obj_num : num_objs) {

        make_obj_data(obj_num, car_nodes, objects_vector);

        vector<int> k_vals = {30};
//        vector<int> k_vals = {1,10,20,30,40};
        for (int k: k_vals) {

            int para1[100] = {0};
            int para2[100] = {0};
            int para3[100] = {0};
            int index = 0;
            tradeoff_level = 0;
            for (cut_level = 0; cut_level <= level_num; cut_level++) {
                para1[index] = tradeoff_level;
                para2[index] = cut_level;
                para3[index++] = 2;
            }
            tradeoff_level = level_num;
            for (cut_level = level_num; cut_level >= 0; cut_level--) {
                para1[index] = tradeoff_level;
                para2[index] = cut_level;
                para3[index++] = 3;
            }
            cout << "obj num: " << obj_num << endl;
            vector<simulate_times> cgraph_simulate_times;
            string method_name = "ours";

            string simulate_file_name = make_out_file_name(method_name, k, obj_num / test_n, "simulate");
            FILE *simulate_file = fopen((input_parameters.output_data_dir + simulate_file_name).c_str(), "w");

            if (simulate_file == NULL) {
                cout << "simulate file is NULL in simulate_toain_query_first()" << endl;
                stop_here();
            }

            // there are a number of (at most $index) algorithms, let's test them one by one
            for (int test_case = 0; test_case < index; test_case++) {
                tradeoff_level = para1[test_case];
                cut_level = para2[test_case];
                mode = para3[test_case];
                simulate_times current_simulate_time = collect_simulate_uber_one_instance(car_nodes, k, mems, obj_num,
                                                                                          num_of_periods,
                                                                                          min_avg_queries_per_period, 3);
                cgraph_simulate_times.push_back(current_simulate_time);
            }

            for (int i = 0; i < cgraph_simulate_times.size(); i++) {
                fprintf(simulate_file, "Combination %d\n", i);
                fprintf(simulate_file, "Query costs\n");
                for (int j = 0; j < cgraph_simulate_times[i].query_times.size(); j++) {
                    for (int z = 0; z < cgraph_simulate_times[i].query_times.size(); z++) {
                        fprintf(simulate_file, "%d \t %d\n", j, cgraph_simulate_times[i].query_times[j][z]);
                    }
                }
                fprintf(simulate_file, "Update costs\n");

                for (int j = 0; j < cgraph_simulate_times[i].update_times.size(); j++) {

                    fprintf(simulate_file, "%d\n", cgraph_simulate_times[i].update_times[j]);

                }
            }
            fclose(simulate_file);

        }

    }

    delete_mems(mems);
    delete[] car_nodes;
}


void simulate_toain_fifo() {
    // Prepare
    test_parameters.arrival_model = 0;
    srand(time(NULL));
    read_road_network();
    set_base_cell_num();
    int max_taxi_num = test_n;
    read_scob_index();
    mem_struct mems;
    allocate_mem(mems, test_n + 1 + max_taxi_num);
    int k_max = 50;
    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;

    int *car_nodes = new int[test_n + 1];
    vector<double> num_objs;

    if (test_parameters.data == 0)
        for (double nums = test_n * 0.001; nums <= test_n * 0.001; nums *= 10)
            num_objs.push_back(nums);
    else
        num_objs.push_back(0);

//    int duration_secs = 30;
    int num_of_updates = 2000;
    int num_of_queries_per_update = 1;
    int simulation_duration=30;
    for (double obj_num : num_objs) {

        make_obj_data(obj_num, car_nodes, objects_vector);

        vector<int> k_vals = {30};

        for (int k: k_vals) {

            int para1[100] = {0};
            int para2[100] = {0};
            int para3[100] = {0};
            int index = 0;
            tradeoff_level = 0;
            for (cut_level = 0; cut_level <= level_num; cut_level++) {
                para1[index] = tradeoff_level;
                para2[index] = cut_level;
                para3[index++] = 2;
            }
            tradeoff_level = level_num;
            for (cut_level = level_num; cut_level >= 0; cut_level--) {
                para1[index] = tradeoff_level;
                para2[index] = cut_level;
                para3[index++] = 3;
            }
            cout << "obj num: " << obj_num << endl;
            vector<simulate_times> cgraph_simulate_times;
            string method_name = "ours";

            string simulate_file_name = make_out_file_name(method_name, k, obj_num / test_n, "simulate");
            FILE *simulate_file = fopen((input_parameters.output_data_dir + simulate_file_name).c_str(), "w");

            if (simulate_file == NULL) {
                cout << "simulate file is NULL in simulate_toain_fifo()" << endl;
                stop_here();
            }

            // there are a number of (at most $index) algorithms, let's test them one by one
            for (int test_case = 0; test_case < index; test_case++) {
                tradeoff_level = para1[test_case];
                cut_level = para2[test_case];
                mode = para3[test_case];
                simulate_times current_simulate_time = collect_simulate_poke_one_instance(k_max, k, car_nodes, mems,
                                                                                          obj_num,
                                                                                          num_of_updates,
                                                                                          num_of_queries_per_update,
                                                                                          simulation_duration);
                cgraph_simulate_times.push_back(current_simulate_time);
            }

            for (int i = 0; i < cgraph_simulate_times.size(); i++) {
                fprintf(simulate_file, "Combination %d\n", i);
                fprintf(simulate_file, "Query costs\n");
                for (int j = 0; j < cgraph_simulate_times[i].query_times.size(); j++) {
                    for(int z =0;z<cgraph_simulate_times[i].query_times[j].size();z++)
                        fprintf(simulate_file, "%d \t %d\n", j, cgraph_simulate_times[i].query_times[j][z]);
                }
                fprintf(simulate_file, "Update costs\n");
                for (int j = 0; j < cgraph_simulate_times[i].update_times.size(); j++) {
                    fprintf(simulate_file, "%d\n", cgraph_simulate_times[i].update_times[j]);
                }
            }
            fclose(simulate_file);

        }

    }

    delete_mems(mems);
    delete[] car_nodes;
}


void simulate_dijk() {
    // Prepare
    test_parameters.arrival_model = 0;
    srand(time(NULL));
    read_road_network();
    int max_taxi_num = test_n;
    mem_struct mems;
    allocate_mem(mems, test_n + 1 + max_taxi_num);

    read_scob_index();
    int *car_nodes = new int[test_n + 1];
    vector<double> num_objs;

    if (test_parameters.data == 0)
        for (double nums = test_n * 0.001; nums <= test_n * 0.001; nums *= 10)
            num_objs.push_back(nums);
    else
        num_objs.push_back(0);

    int duration_secs = 10;
    for (double obj_num : num_objs) {

        make_obj_data(obj_num, car_nodes, objects_vector);

        vector<int> k_vals = {1, 10, 20, 30, 40};

        for (int k: k_vals) {


            cout << "obj num: " << obj_num << endl;
            vector<long> dijk_simulate_times;
            string method_name = "dijk";

            string simulate_file_name = make_out_file_name(method_name, k, obj_num / test_n, "simulate");
            FILE *simulate_file = fopen((input_parameters.output_data_dir + simulate_file_name).c_str(), "w");

            if (simulate_file == NULL) {
                cout << "simulate file is NULL in simulate_dijk()" << endl;
                stop_here();
            }
            long duration = 0;

            while (true) {
                int query=get_valid_random_node();
                long start = clock();
                naiveKNN(k, objects_vector, query, mems.dist, mems.visited, mems.q, car_nodes);
                duration += clock() - start;
                dijk_simulate_times.push_back(clock() - start);
                if (duration > duration_secs*MICROSEC_PER_SEC) break;
            }

            for (int j = 0; j < dijk_simulate_times.size(); j++) {
                fprintf(simulate_file, "%ld\n", dijk_simulate_times[j]);
            }

            fclose(simulate_file);

        }

    }

    delete_mems(mems);
    delete[] car_nodes;
}

#endif //SOB_SIMULATION_H
