/*
 * TestMethods.h
 *
 *  Created on: 2016年8月29日
 *      Author: lostrong
 */

#ifndef SRC_TESTMETHODS_H_
#define SRC_TESTMETHODS_H_

#include<vector>
//#include "PivotKNN.h"
#include "TestToain.h"
#include "../scob/ScobQuery.h"
#include "../util/TimeStat.h"
#include "../RoadKNNUtilities.h"
#include <iomanip>
#include "../util/MemoryAllocation.h"

using namespace std;

vector<int> get_rand_cars(int num) {
    vector<int> results;
    srand(time(NULL));
    int *used = new int[test_n + 1];
    memset(used, 0, sizeof(int) * (test_n + 1));
    for (int i = 0; i < num; ++i) {
        int node;
        do {
            node = get_rand_node();
        } while (used[node]);
        results.push_back(node);
        used[node] = 1;
    }
    delete[] used;
    return results;
}


void update_cars_naive(int k_max, vector<CarPos> &carposes, vector<CarPos> &allcarposes, int update_item) {
    for (int i = 0; i < carposes.size(); i++)
        if (carposes[i].id == allcarposes[update_item].id) {
            carposes[i] = allcarposes[update_item];
//			break;
        }
}

void update_cars_hier_node(int k_max, int k, vector<int> &node_ran_objs, int update_item, mem_struct &mems) {
    int n1 = node_ran_objs[update_item];


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

    remove_kdnn_object(n1, k_max, k, mems, rev_shortcuts_arr, 10);

    update_kDNN_object(0, n1, k_max, k, mems.dist, mems.visited, mems.q);

}

void update_cars_hier(int k_max, int k, unordered_map<int, vector<PointCar> > &pointcar_hash,
                      unordered_map<int, CarPos> &carposes_map, vector<CarPos> &allcarposes, int update_item,
                      DijkstraQueue *q1,
                      long *dis1, int *visited1) {
//	cout << update_item << endl;
    int n1 = allcarposes[update_item].node1;
    int n2 = allcarposes[update_item].node2;
    int d1 = allcarposes[update_item].dis1;
    int d2 = allcarposes[update_item].dis2;
    int id = allcarposes[update_item].id;
    int ori_n1 = carposes_map[id].node1;
    int ori_n2 = carposes_map[id].node2;
    vector<PointCar> ori_v1 = pointcar_hash[ori_n1];
    vector<PointCar> ori_v2 = pointcar_hash[ori_n2];
    for (vector<PointCar>::iterator it = ori_v1.begin(); it != ori_v1.end(); it++) {

        if (it->carid == id) {
            it = ori_v1.erase(it);
            vector<KNode> list = hier_local_knn_arr[ori_n1];
            for (vector<KNode>::iterator iter = list.begin(); iter != list.end(); iter++) {
                if (iter->id == id + test_n) {
                    iter = list.erase(iter);
                    break;
                }
            }

            hier_local_knn_arr[ori_n1] = list;
            break;
        }

    }
    pointcar_hash[ori_n1] = ori_v1;

    int cnt2 = 0;
    for (vector<PointCar>::iterator it = ori_v2.begin(); it != ori_v2.end(); cnt2++, it++) {

        if (it->carid == id) {
            it = ori_v2.erase(it);
            vector<KNode> list = hier_local_knn_arr[ori_n2];
            for (vector<KNode>::iterator iter = list.begin(); iter != list.end(); iter++) {
                if (iter->id == id + test_n) {
                    iter = list.erase(iter);

                    break;
                }
            }
            hier_local_knn_arr[ori_n2] = list;
            break;
        }
    }
    pointcar_hash[ori_n2] = ori_v2;
    remove_hier_local_knn_carpos(carposes_map[id], k_max, k, pointcar_hash, q1, dis1, visited1);

    vector<PointCar> v1 = pointcar_hash[n1];
    vector<PointCar> v2 = pointcar_hash[n2];

    int i = 0;
    for (i = 0; i < v1.size(); i++) {
        if (d1 < v1[i].dis) {
            break;
        }
    };
    v1.push_back(PointCar(id, d1));
    int size = v1.size();

    for (int j = size - 1; j > i; j--)
        v1[j] = v1[j - 1];

    v1[i] = PointCar(id, d1);
    pointcar_hash[n1] = v1;

    i = 0;
    for (i = 0; i < v2.size(); i++) {
        if (d2 < v2[i].dis) {
            break;
        }
    }

    v2.push_back(PointCar(id, d2));
    size = v2.size();
    for (int j = size - 1; j > i; j--)
        v2[j] = v2[j - 1];
    PointCar tmpnode2(id, d2);
    v2[i] = tmpnode2;

    pointcar_hash[n2] = v2;

    carposes_map[id] = allcarposes[update_item];
//	cout<<"here4"<<endl;
    update_hier_local_knn_carpos_edgeversion(allcarposes[update_item], k_max, dis1, visited1, q1);

}


void verify_two_results(vector<KNode> &r1, vector<KNode> &r2, int update, int query, int k, mem_struct &mems) {

    if (r1.size() != r2.size()) {
        cout << "update: " << update << endl;
        cout << "r1 size != r2 size..." << endl;
        cout << query << " " << r1.size() << " " << r2.size() << endl;
        stop_here();
    }
    for (int i = 0; i < k && i < r1.size(); ++i) {
        if (r1[i].dis != r2[i].dis) {
            cout << "update: " << update << endl;
            int id1, id2;
            if (test_parameters.version == NODE_VERSION_RAN) {
                id1 = r1[i].id;
                id2 = r2[i].id;
            }
            if (test_parameters.version == EDGE_VERSION_RAN) {
                id1 = carposes_map[r1[i].id - test_n].node1;
                id2 = carposes_map[r2[i].id - test_n].node1;
            }
            //					cout << "carpose to be updated: " << allcarposes[update].id
            //							<< " " << allcarposes[update].node1 << " "
            //							<< allcarposes[update].dis1 << " "
            //							<< allcarposes[update].node2 << " "
            //							<< allcarposes[update].dis2 << endl;
            cout << query << " " << id1 << " " << i << " not equal!!" << endl;
            cout << get_distance_st(query, id1) << endl;
            cout << TestShortestPath(query, id1, node_order, shortcuts, mems, INT_MAX) << endl;
            cout << hierSP_rev(query, id1, node_order, shortcuts, mems, INT_MAX) << endl;
            cout << test_shortest_path_bidirection(query, id1, node_order, shortcuts, mems, INT_MAX) << endl;
            for (int j = 0; j < r1.size(); ++j) {
                if (test_parameters.version == NODE_VERSION_RAN) {
                    cout << r1[j].id << " " << r1[j].dis

                         << endl;
                    cout << r2[j].id << " " << r2[j].dis << endl << endl;
                }
                if (test_parameters.version == EDGE_VERSION_RAN) {
                    cout << r1[j].id - test_n << " " << r1[j].dis

                         << endl;
                    cout << r2[j].id - test_n << " " << r2[j].dis << endl << endl;
                }
            }

            stop_here();
        }

    }
}

void test_toain_free_one_instance(int k_max, int k, vector<times> &cgraph_times, vector<times> &dijk_times,
                                  vector<int> &queries, int *car_nodes, mem_struct &mems, double obj_num) {
    times cgraph_current_times;
    times dijk_current_times;
    DijkstraQueue *q = mems.q;
    long *dist = mems.dist;
    int *visited = mems.visited;

    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;
    //	init_hier_local_knn_carpos(carposes, k_max);


    vector<long> dijk_each_query_time;
    vector<long> cgraph_each_query_time;
    vector<long> cgraph_each_update_time;
    /*
     * cal init time for C-Graph, it is also the time spent on updating all the objects
     */
    long start_init = clock();
    for (int i = 0; i <= test_n; i++) {
        hier_local_knn_arr[i].clear();
    }
    if (test_parameters.arrival_model == 1) {
        cgraph_each_update_time = update_kdnn_list(1, k_max, 0, test_parameters.version, q, dist, visited,
                                                   objects_vector);
    } else {
        update_kdnn_list(1, k_max, 0, test_parameters.version, q, dist, visited, objects_vector);
    }
    long cgraph_init_time = clock() - start_init;
    cout << "cgraph init time: " << cgraph_init_time << endl;
    //			dijk_init_per_obj_times.push_back(cgraph_init_time/(test_n*ratio));
    cout << "cgraph init time (per object): " << cgraph_init_time / obj_num << endl;
    cgraph_current_times.init = cgraph_init_time;
    //	cgraph_init_per_obj_times.push_back(cgraph_init_time / (test_n * ratio));

    long start_init_dijk = clock();
    update_dijkstra(edge_car_hash);
    long dijk_init_time = clock() - start_init_dijk;
    cout << "dijkstra init time: " << dijk_init_time << endl;
    dijk_current_times.init = dijk_init_time;
    long update1 = 0;
    int num = queries.size();
    // for each setting, we test "num" times, num by default = 200


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

    for (int update = 0; update < num; update++) {
        int query = queries[update];

        /*
         * Method: Dijkstra; Purpose: estimate the update time
         */
        //			long start1 = clock();
        //			update_cars_naive(k_max, carposes, allcarposes, update);
        //			update1 += clock() - start1;
        //
        /*
         * Method: C-Graph; Purpose: estimate the update time for individual update
         */
        int action = 0;
        int affected_node;
        if (test_parameters.arrival_model == 0) {            // pokemon-go
            long start2 = clock();

            if (test_parameters.version == NODE_VERSION_RAN) {
//				int update_item = get_rand_node() % (objects_vector.size());
//				update_cars_hier_node(k_max, objects_vector, update_item, q, dist, visited, marked_first);

                int pnum = rand()%100;
                if (pnum<50) {
                    remove_kdnn_object(objects_vector[0], k_max, k, mems, rev_shortcuts_arr, 10);
                    action = -1;            //remove
                    affected_node = objects_vector[0];

                } else {
                    int inserted_node_id;
                    do {
                        inserted_node_id = get_rand_node();
                    } while (car_nodes[inserted_node_id]);
                    update_kDNN_object(0, inserted_node_id, k_max, k, dist, visited, q);
                    action = 1;
                    affected_node = inserted_node_id;

                }

//				cout<<"update finished"<<endl;
            }
            if (test_parameters.version == EDGE_VERSION_RAN) {
                update_cars_hier(k_max, k, pointcar_hash, carposes_map, allcarposes, update, q, dist, visited);
            }
            cgraph_each_update_time.push_back(clock() - start2);
        }


        if (action == 1) {
            objects_vector.push_back(affected_node);
            car_nodes[affected_node] = 1;
        }
        if (action == -1) {
            objects_vector.erase(objects_vector.begin());
//				cout<<"objects_vector size: "<<objects_vector.size()<<endl;
            car_nodes[affected_node] = 0;
        }

        /*
         * Method: Dijkstra; Purpose: estimate the query time
         */
        long start = clock();
        vector<KNode> r1;
        if (test_parameters.version == NODE_VERSION_RAN) {
            r1 = naiveKNN(k, objects_vector, query, dist, visited, q, car_nodes);
        }
        if (test_parameters.version == EDGE_VERSION_RAN) {
            r1 = naiveKNN_carpos(k, carposes, edge_car_hash, query, dist, visited, q);
        }
        dijk_each_query_time.push_back(clock() - start);

        /*
         * Method: C-Graph; Purpose: estimate the query time
         */
        start = clock();
        vector<KNode> r2;
        if (test_parameters.version == NODE_VERSION_RAN) {
            r2 = scob_KNN_query(k, hier_local_knn_arr, query, dist, visited, q);
        }

        if (test_parameters.version == EDGE_VERSION_RAN) {
            r2 = scob_KNN_query_edge(k, hier_local_knn_arr, query, dist, visited, q);
        }
        cgraph_each_query_time.push_back(clock() - start);
        /*
         * Verify whether the results are correct
         */
        if (test_parameters.verify) {
            verify_two_results(r1, r2, update, query, k, mems);
        }

    }
    dijk_current_times.query = get_mean(dijk_each_query_time);
    dijk_current_times.query_variance = get_var(dijk_each_query_time);
    cgraph_current_times.query_variance = get_var(cgraph_each_query_time);
    cgraph_current_times.update_variance = get_var(cgraph_each_update_time);
    cgraph_current_times.query = get_mean(cgraph_each_query_time);
    cgraph_current_times.update = get_mean(cgraph_each_update_time);

    //		estimate_hito(file_query, cgraph_each_query_time, 20);
    //		estimate_hito(file_update, cgraph_each_update_time, 20);
    //		for (long qtime : cgraph_each_query_time)
    //			fprintf(file_query, "%ld\n", qtime);
    //		fprintf(file_query, "\n\n");
    //
    //		for (long qtime : cgraph_each_update_time)
    //			fprintf(file_update, "%ld\n", qtime);
    //		fprintf(file_update, "\n\n");

    dijk_times.push_back(dijk_current_times);
    cgraph_times.push_back(cgraph_current_times);
    cout << "cgraph query time: " << cgraph_current_times.query << endl;
    //			//note: num times of calling the udpate methods, equals to estimating num/2 times removal, num/2 times
    //			//insertions and num/2 times updates, in total 1.5num times.
    //			//		cgraph_update_per_obj_times.push_back(cgraph_total_update_time / (num * 1.5));
    cout << "cgraph update time: " << cgraph_current_times.update << endl;
    			cout << "dijkstra query time: " << dijk_current_times.query << endl;



}


void test_toain_free() {
    // Prepare
    srand(time(NULL));
    read_road_network();
    set_base_cell_num();
    int max_taxi_num = test_n;
    int k_max = 50;
    read_scob_index();
    mem_struct mems;
    allocate_mem(mems, test_n + 1 + max_taxi_num);
    vector<int> accumulated_comb_ids;
    level_num = get_level_num_with_dataspace();
    cout << "level_num: " << level_num << endl;

    // Read queries
    int num = 200;

    vector<int> queries_from_file = read_queries(network_name.c_str(), num);
    cout << "read query ends..." << endl;
    long tuning_time = clock();
    vector<int> queries;
    for (int eta = 0; eta < num; eta++) {
        queries.push_back(queries_from_file[eta]);
    }
    int *car_nodes = new int[test_n + 1];
    vector<double> num_objs;


    for (double nums = test_n * 0.001; nums <= test_n * 0.1; nums *= 10)
        num_objs.push_back(nums);

//    if (test_parameters.arrival_model == 1)
//        file = fopen("NW_POI_Ours_k1.txt", "w");
//    else
//        file = fopen("NW_POI_Ours_k1_poke++.txt", "w");

    for (double obj_num : num_objs) {
        if (test_parameters.version == NODE_VERSION_RAN) {
            obj_num = make_obj_data(obj_num, car_nodes, objects_vector);
            cout << (int) (obj_num) << " objects" << endl;
        } else {
            carposes = carposes_random;

            if (test_parameters.version == EDGE_VERSION_RAN)
                // make carposes_random
                read_rand_car_positions((int) (obj_num), carposes);
            if (test_parameters.version == EDGE_VERSION_REAL)
                read_car_positions("init_car_positions.txt", carposes, carposes_map, pointcar_hash);
        }
//        vector<int> k_vals = {1, 10, 20, 30, 40};
        vector<int> k_vals = {1};
        for (int k:k_vals) {

            string file_name = make_out_file_name("ours", k, obj_num, "analytical");
            FILE *file = fopen((input_parameters.output_data_dir + file_name).c_str(), "w");
            int para1[100] = {0};
            int para2[100] = {0};
            int para3[100] = {0};
            int index = 0;
//            tradeoff_level = 0;
//            for (cut_level = 0; cut_level <= level_num; cut_level++) {
//                para1[index] = tradeoff_level;
//                para2[index] = cut_level;
//                para3[index++] = 2;
//            }
            tradeoff_level = level_num;
            for (cut_level = level_num; cut_level >= 0; cut_level--) {
                para1[index] = tradeoff_level;
                para2[index] = cut_level;
                para3[index++] = 3;
            }
            cout << "obj num: " << obj_num << endl;
            vector<times> cgraph_times;
            vector<times> dijk_times;


            // there are a number of (at most $index) algorithms, let's test them one by one
            for (int test_case = 0; test_case < index; test_case++) {
                tradeoff_level = para1[test_case];
                cut_level = para2[test_case];
                mode = para3[test_case];
                test_toain_free_one_instance(k_max, k, cgraph_times, dijk_times, queries, car_nodes, mems, obj_num);

            }

            cout << "Tuning time: " << clock() - tuning_time << endl;
            cout << "detailed data:" << endl;

            cout << setw(12) << "Combination: " << setw(2) << setw(14) << "query" << setw(12) << "variance" << setw(12)
                 << "update" << setw(16) << "variance" << endl;
            for (int i = 0; i < dijk_times.size(); i++) {
                cout << setw(10) << "Combination: " << setw(2) << i << setw(12) << cgraph_times[i].query << setw(12)
                     << cgraph_times[i].query_variance << setw(12) << cgraph_times[i].update << setw(16)
                     << cgraph_times[i].update_variance << endl;
            }

//		stop_here();

            if (test_parameters.arrival_model == 1)
                analytical_query_first_all_parameters(cgraph_times, dijk_times, obj_num, file,
                                                      accumulated_comb_ids);
            else
                analytical_fifo_all_parameters(cgraph_times, dijk_times, obj_num, file, accumulated_comb_ids);



            fclose(file);
        }
        if(test_parameters.data!=0)
            break;

    }

    sort(accumulated_comb_ids.begin(), accumulated_comb_ids.end());
    int cnt = 0;
    int pre = -1;
    cout << endl << "id  cnt:" << endl;
    for (int id : accumulated_comb_ids) {
        if (id != pre) {
            if (pre != -1) {
                cout << pre << " " << cnt << endl;
            }
            pre = id;
            cnt = 1;

        } else {
            cnt++;
        }
    }
    if (cnt > 0) {
        cout << pre << " " << cnt << endl;
    }
    delete_mems(mems);
    delete[] car_nodes;
}

#endif /* SRC_TESTMETHODS_H_ */
