//
// Created by lostrong on 17/3/8.
//

#ifndef SOB_TRAINING_H
#define SOB_TRAINING_H

#include "../para/TestParameter.h"
#include "../para/GlobalVariables.h"
#include "../RoadKNNUtilities.h"
#include "../scob/ScobIndex.h"
#include "../util/DataFunctions.h"
#include "../util/TimeStat.h"
#include "../scob/ScobQuery.h"

/**
 *
 * @param k
 * @param queries
 * @param mems
 * @param obj_num
 * @param comb_id
 * @return
 * The training process is simple. We generate num ( = 500) queries, and estimate the mean query/update time etc. Then
 * use formulae to compute the throughput.
 */
times collect_training_query_first_one_instance(int k,
                                                vector<int> &queries, mem_struct &mems, double obj_num, int comb_id) {
    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;
    times cgraph_current_times;
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;
    vector<long> cgraph_each_query_time;
    vector<long> cgraph_each_update_time;
//    long duration;
    long start_init = clock();
    for (int i = 0; i <= test_n; i++) {
        hier_local_knn_arr[i].clear();
    }
    int training_periods = 1;
    long cgraph_init_time = 0;
    for(int i = 0; i< training_periods; i++) {
        int is_init = 1;
        cout<<"starting init..."<<endl;
        // when is_init = 1, we initialize the kDNN list, and k_max is set to k, k is set to 0 (a not important value);
        cgraph_each_update_time = update_kdnn_list(is_init, k, 0, test_parameters.version, q, dist, visited,
                                                   objects_vector);
        cgraph_init_time += clock() - start_init;
        cout << "cgraph init time: " << cgraph_init_time << endl;
        cout << "cgraph init time (per object): " << (clock() - start_init) / obj_num << endl;


        int num = queries.size();
        // for each setting, we test "num" times, num by default = 500
        for (int update = 0; update < num; update++) {
            int query = queries[update];
            long start = clock();
            vector<KNode> r;
            r = scob_KNN_query(k, hier_local_knn_arr, query, dist, visited, q);
            cgraph_each_query_time.push_back(clock() - start);
        }
    }
    cgraph_current_times.init = cgraph_init_time/training_periods;
    cgraph_current_times.k = k;
    cgraph_current_times.network_name = network_name;
    if (test_parameters.data == DATA_RANDOM)
        cgraph_current_times.obj_name = "ran";
    else
        cgraph_current_times.obj_name = "real";
    cgraph_current_times.method_name = "ours";
    cgraph_current_times.obj_ratio = process(obj_num * 1.0 / test_n, 6);
    cgraph_current_times.query_variance = get_var(cgraph_each_query_time);
    cgraph_current_times.update_variance = get_var(cgraph_each_update_time);
    cgraph_current_times.query = get_mean(cgraph_each_query_time);
    cgraph_current_times.update = get_mean(cgraph_each_update_time);
    cgraph_current_times.is_valid = 1;
    cgraph_current_times.scob_configuration_id = comb_id;
    cout << "cgraph query time: " << cgraph_current_times.query << endl;
    cout << "cgraph update time: " << cgraph_current_times.update << endl;
    return cgraph_current_times;
}


// Fixing k_max = 50
times collect_training_fifo_one_instance(int k_max, int k,
                                         vector<int> &queries, int *&object_nodes, mem_struct &mems, double obj_num) {
    int queries_per_update = 50;
    for (int u = 0; u < test_n; u++)
        rev_shortcuts_arr[u].clear();

    for (int u = 0; u < test_n; u++) {

        vector<KNode> *nei;
        if (mode == QUERY_FAVOR_MODE) {
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
    times toain_current_times;
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;
    // for certain configurations, update can be very costly, and these configurations are not useful. So we
    // set a limit on the training time to filter out these configurations.
    int training_duration_limits = 60 * MICROSEC_PER_SEC;
    long duration = 0;
    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;
    vector<long> toain_each_query_time;
    vector<long> toain_each_update_time;
    vector<int> non_object_nodes;
    long start_init = clock();
    for (int i = 0; i <= test_n; i++) {
        hier_local_knn_arr[i].clear();
    }

    update_kdnn_list(1, k_max, 0, test_parameters.version, q, dist, visited, objects_vector);
    for (int i = 1; i < test_n; i++) {
        if (object_nodes[i] == 0) {
            non_object_nodes.push_back(i);
        }
    }

//    for (int i = 0; i < SHUFFLE_NUM; i++) {
//        // shuffle non objects
//        int z1 = rand() % non_object_nodes.size();
//        int z2 = rand() % non_object_nodes.size();
//        int tmp = non_object_nodes[z1];
//        non_object_nodes[z1] = non_object_nodes[z2];
//        non_object_nodes[z2] = tmp;
//
//        // shuffle objects
//        z1 = rand() % objects_vector.size();
//        z2 = rand() % objects_vector.size();
//        tmp = objects_vector[z1];
//        objects_vector[z1] = objects_vector[z2];
//        objects_vector[z2] = tmp;
//    }
    shuffle(objects_vector);
    shuffle(non_object_nodes);
    long toain_init_time = clock() - start_init;
    cout << "cgraph init time: " << toain_init_time << endl;
    cout << "cgraph init time (per object): " << toain_init_time / obj_num << endl;
    toain_current_times.init = toain_init_time;


//    int num = queries.size();
    int num = 10000;
    // for each setting, we test "num" times, num by default = 200
    for (int update = 0; update < num; update++) {
//        int query = queries[update];
        int pnum = rand() % 100;
        if ((pnum < 50 && (objects_vector.size() > 0)) || (non_object_nodes.size() == 0)) {
            // since objects_vector has been shuffled, selecting the first node is a random selection
            int removeNode = objects_vector[0];
            long start1 = clock();
            int each_removal_duration_limit = 30;// set each removal not longer than 30 seconds
            int succeed_within_limit = remove_kdnn_object(removeNode, k_max, k, mems, rev_shortcuts_arr,
                                                          each_removal_duration_limit);
            if (succeed_within_limit == 0) break;
            toain_each_update_time.push_back(clock() - start1);
            duration += clock() - start1;
            objects_vector.erase(objects_vector.begin());
            object_nodes[removeNode] = 0;
            non_object_nodes.push_back(removeNode);

        } else {
            long start2 = clock();
            int inserted_node_id = non_object_nodes[0];
            int is_init = 0;// not init, update an object
            update_kDNN_object(is_init, inserted_node_id, k_max, k, dist, visited, q);
            toain_each_update_time.push_back(clock() - start2);
            duration += clock() - start2;
            non_object_nodes.erase(non_object_nodes.begin());
            objects_vector.push_back(inserted_node_id);
            object_nodes[inserted_node_id] = 1;
        }


        vector<int> query_nodes;
        for (int zz = 0; zz < queries_per_update; zz++) {
            int query;
            do {
                query = rand() % test_n;
            } while (valid_node_hash[query] == 0);
            query_nodes.push_back(query);
        }
//        vector<KNode> r2;
        for (int zz = 0; zz < queries_per_update; zz++) {
//            cout<<"query: "<<query_nodes[zz]<<endl;

//            vector<KNode> r1 = naiveKNN(k, objects_vector, query_nodes[zz], dist, visited, q, object_nodes);
            long start = clock();
//            vector<KNode> r2=
            scob_KNN_query(k, hier_local_knn_arr, query_nodes[zz], dist, visited, q);
            toain_each_query_time.push_back(clock() - start);
            duration += clock() - start;
//            verify_two_results(r1, r2, update, query_nodes[zz], k, mems);
        }


        if (duration > training_duration_limits)
            break;
    }

    toain_current_times.
            k = k;
    if (test_parameters.data == DATA_RANDOM)
        toain_current_times.
                obj_name = "ran";
    else
        toain_current_times.
                obj_name = "real";
    toain_current_times.
            method_name = "ours";
    toain_current_times.
            obj_ratio = process(obj_num * 1.0 / test_n, 6);
    toain_current_times.
            query_variance = get_var(toain_each_query_time);
    toain_current_times.
            update_variance = get_var(toain_each_update_time);
    toain_current_times.
            query = get_mean(toain_each_query_time);
    toain_current_times.
            update = get_mean(toain_each_update_time);
    return
            toain_current_times;

}

void train_toain_fifo() {
    // Prepare
    test_parameters.arrival_model = QUEUE_POLICY_FIFO;
    srand(time(NULL));
    read_road_network();
    set_base_cell_num();
    int max_taxi_num = test_n;
    int k_max = 50;// for fifo model, we always preset max_k = 50 for the testing on paper;
    read_scob_index();
    mem_struct mems;
    allocate_mem(mems, test_n + 1 + max_taxi_num);

    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;

    // Read queries
    int num = 1000;
    vector<int> queries;
    for (int eta = 0; eta < num; eta++) {
        int query;
        do {
            query = get_rand_node();
        } while (valid_node_hash[query] == 0);
        queries.push_back(query);
    }
    int *car_nodes = new int[test_n + 1];
    vector<double> num_objs;


    if (test_parameters.data == DATA_RANDOM) {
        vector<double> ratios = {DEFAULT_RATIO};
//        vector<double> ratios = {5000.0/test_n};
        for (double ratio: ratios) {
            num_objs.push_back(test_n * ratio);
            double num = test_n*ratio;
            cout<<"here: "<<num/test_n<<endl;
        }
    } else {
        // If use real objects, then randomly push a number into the object set; pushing in what number is not important
        num_objs.push_back(1);
    }

    vector<int> k_vals = {1, 10, 20, 30, 40};
//    vector<int> k_vals = {DEFAULT_K};
//    vector<int> k_vals = {9};// case study 1
    cout<<"network_name: "<<network_name<<endl;
    for (double obj_num : num_objs) {
        int obj_num_2 = make_obj_data(obj_num, car_nodes, objects_vector);
        if(test_parameters.data != DATA_RANDOM) {
            obj_num=obj_num_2;
        }
        for (int k: k_vals) {
            vector<ScobConfiguration> configurations;
            if(network_name.compare("BJ")==0)
                configurations = ScobConfiguration::generate_all_update_favor_configurations(
                    level_num);
            else
                configurations = ScobConfiguration::generate_all_configurations(
                        level_num);

            cout << "obj num: " << obj_num << endl;

            vector<times> toain_train_times;
            string method_name = "ours";
            cout<<"obj ratio: "<<obj_num/test_n<<endl;
            string tune_file_name = make_out_file_name(method_name, k, obj_num / test_n, "tune");
            FILE *tune_file = fopen((input_parameters.output_data_dir + tune_file_name).c_str(), "w");
            test_valid_file(tune_file, (input_parameters.output_data_dir + tune_file_name).c_str(), "train_toain_fifo");
            string train_file_name = make_out_file_name(method_name, k, obj_num / test_n, "training");
            FILE *train_file = fopen((input_parameters.output_data_dir + train_file_name).c_str(), "w");
            test_valid_file(train_file, (input_parameters.output_data_dir + train_file_name).c_str(),
                            "train_toain_fifo");
            long tuning_time = clock();
            // there are a number of (at most $index) configurations, let's test them one by one
            for (int test_case = 0; test_case < configurations.size(); test_case++) {
                tradeoff_level = configurations[test_case].tradeoff_level;
                cut_level = configurations[test_case].cut_level;
                mode = configurations[test_case].mode;
                times cgraph_current_time = collect_training_fifo_one_instance(k_max, k, queries, car_nodes, mems,
                                                                               obj_num);
                toain_train_times.push_back(cgraph_current_time);
                cout << "finish case: " << test_case << endl;

            }
            cout << "finish all test cases" << endl;
            fprintf(tune_file, "%d \t %lf \t %ld\n", k, obj_num, clock() - tuning_time);
            fprintf(train_file, "k \t obj_ratio \t Combination \t query \t variance \t update \t variance\n");

            for (int i = 0; i < toain_train_times.size(); i++) {
                fprintf(train_file, "%d \t %lf \t %d \t %lf \t %lf \t %lf \t %lf \n", toain_train_times[i].k,
                        toain_train_times[i].obj_ratio, i, toain_train_times[i].query,
                        toain_train_times[i].query_variance,
                        toain_train_times[i].update, toain_train_times[i].update_variance);
            }
            fclose(train_file);
            fclose(tune_file);

        }
    }
    cout << "start release memory..." << endl;
    delete_mems(mems);
    delete[] car_nodes;
}


void train_toain_query_first() {
    // Prepare
    test_parameters.arrival_model = QUEUE_POLICY_QUERY_FIRST;
    srand(time(NULL));
    read_road_network();
    set_base_cell_num();
    int max_taxi_num = test_n;
    read_scob_index();
    mem_struct mems;
    allocate_mem(mems, test_n + 1 + max_taxi_num);

    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;

    vector<int> k_vals = {DEFAULT_K};
    vector<double> obj_ratios;
    // Make queries
    int num = 1000;
    int query;
    vector<int> queries;
    for (int i = 0; i < num; i++) {
        do {
            query = get_rand_node();
            queries.push_back(query);
        } while (valid_node_hash[query] == 0);
    }
//    vector<int> queries= read_queries(network_name.c_str(), num);
    cout << "read query ends..." << endl;
    int *car_nodes = new int[test_n + 1];

    // Make objects
    vector<double> num_objs;
    if (test_parameters.data == DATA_RANDOM)
        obj_ratios = {DEFAULT_RATIO};
    else
        obj_ratios = {1.0};


    for (double ratio : obj_ratios)
        num_objs.push_back(test_n * ratio);

    // for training of query first model, we use k and object ratio
    for (double obj_num : num_objs) {
        obj_num = make_obj_data(obj_num, car_nodes, objects_vector, 4);
        for (int k: k_vals) {
            vector<ScobConfiguration> configurations = ScobConfiguration::generate_all_configurations(level_num);
            cout << "obj num: " << obj_num << endl;
            vector<times> cgraph_train_times;
            string method_name = "ours";
            string tune_file_name = make_out_file_name(method_name, k, obj_num / test_n, "tune");
            FILE *tune_file = fopen((input_parameters.output_data_dir + tune_file_name).c_str(), "w");
            test_valid_file(tune_file, (input_parameters.output_data_dir + tune_file_name).c_str(),
                            "train_toain_query_first");
            string train_file_name = make_out_file_name(method_name, k, obj_num / test_n, "training");
            FILE *train_file = fopen((input_parameters.output_data_dir + train_file_name).c_str(), "w");
            test_valid_file(train_file, (input_parameters.output_data_dir + train_file_name).c_str(),
                            "train_toain_query_first");

            long tuning_time = clock();
            // there are a number of (at most $index) algorithms, let's test them one by one
            for (int test_case = 0; test_case < configurations.size(); test_case++) {
                tradeoff_level = configurations[test_case].tradeoff_level;
                cut_level = configurations[test_case].cut_level;
                mode = configurations[test_case].mode;
                times current_train_time = collect_training_query_first_one_instance(k, queries, mems, obj_num,
                                                                                     test_case);
                cgraph_train_times.push_back(current_train_time);
            }
            fprintf(tune_file, "%d \t %lf \t %ld\n", k, obj_num, clock() - tuning_time);
            fprintf(train_file, "k \t obj_ratio \t Combination \t query \t variance \t update \t variance\n");

            for (int i = 0; i < cgraph_train_times.size(); i++) {
                fprintf(train_file, "%d \t %lf \t %d \t %f \t %lf \t %f \t %lf \n", cgraph_train_times[i].k,
                        cgraph_train_times[i].obj_ratio, i, cgraph_train_times[i].query,
                        cgraph_train_times[i].query_variance,
                        cgraph_train_times[i].update, cgraph_train_times[i].update_variance);
            }
            fclose(train_file);
            fclose(tune_file);

        }

    }


    delete_mems(mems);
    delete[] car_nodes;
}


#endif //SOB_TRAINING_H
