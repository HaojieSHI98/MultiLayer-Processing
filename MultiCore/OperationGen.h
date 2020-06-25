//
// Created by lostrong on 18/7/1.
//

#ifndef TOAIN_OPERATIONGEN_H
#define TOAIN_OPERATIONGEN_H
#include "MultiCoreVariables.h"
#include "../queue/GeneralQueueSimulator.h"
#include "../MultiCore/Simulation.h"
#include "ArrivalData.h"

vector<int> generate_arrival_nodes(vector<std::pair<double, int> >& full_list, int begin, int end,
long &total_queries_plan, long &total_updates_plan){

    total_queries_plan=0;
    total_updates_plan=0;
//    int rand_idx_update=0;
//    int rand_idx_query=0;
    vector<int> object_list;
    vector<int> non_object_list;
    vector<int> arrival_node_list;
    for (int i = begin; i < end; i++) {
//            if(vStart[i])
        non_object_list.push_back(i);
    }

    for(int i = 0;i<full_list.size();i++){
        pair<double, int> event = full_list[i];
        if (event.second == INSERT && object_list.size() >= test_n){
            arrival_node_list.push_back(-1);
        }
        if(event.second == DELETE && object_list.size() <= 0){
            arrival_node_list.push_back(-1);

        }
        if (event.second == INSERT && object_list.size() < test_n) {
            total_updates_plan++;
//                int index = global_random_numbers[rand_idx_update] % non_object_list.size();
            int index = rand() % non_object_list.size();
//            rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;
            int non_object_node = non_object_list[index];
            int last_index = non_object_list.size() - 1;

            int tmp = non_object_list[index];
            non_object_list[index] = non_object_list[last_index];
            non_object_list[last_index] = tmp;
            non_object_list.pop_back();

            object_list.push_back(non_object_node);
            arrival_node_list.push_back(non_object_node);

        }
        if (event.second == DELETE && object_list.size() > 0) {
            total_updates_plan++;
//                int index = global_random_numbers[rand_idx_update] % object_list.size();
            int index = rand() % object_list.size();
//            rand_idx_update=(rand_idx_update+1)%rand_length;
            int object_node = object_list[index];
            int last_index = object_list.size() - 1;
            int tmp = object_list[index];
            object_list[index] = object_list[last_index];
            object_list[last_index] = tmp;
            object_list.pop_back();
            non_object_list.push_back(object_node);
            arrival_node_list.push_back(object_node);

        }
        if (event.second == QUERY) {
            total_queries_plan++;
            int query_node;
            do {
                query_node = rand();
//                rand_idx_query=(rand_idx_query+1)%rand_length;
                query_node = query_node % (test_n-1)+1;

            } while (vStart[query_node] == 0);
            arrival_node_list.push_back(query_node);


        }

    }
    return arrival_node_list;

}

ArrivalData* generate_all_operations(double query_rate, double insert_rate,
                                                  double delete_rate, int simulation_time, int begin_node, int end_node) {
    srand(time(NULL));
    vector <std::pair<double, int>> full_list;
    for (int i = 0; i < multiTestPara.init_objects; i++) {
        full_list.push_back(make_pair(0.0, INSERT));
    }
    vector <std::pair<double, int>> append_list = make_online_query_update_list(query_rate, insert_rate,
                                                                                delete_rate,
                                                                                simulation_time);
    for (pair<double, int> &item : append_list)
        full_list.push_back(item);
    long total_queries_plan=0;
    long total_updates_plan=0;
    vector<int> arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node, total_queries_plan, total_updates_plan);
    cout << "full_list made..." << endl;
    cout << "full list size: " << full_list.size() << endl;
    ArrivalData* ad = new ArrivalData(full_list, arrival_nodes, total_queries_plan, total_updates_plan);
    return ad;
}

#endif //TOAIN_OPERATIONGEN_H
