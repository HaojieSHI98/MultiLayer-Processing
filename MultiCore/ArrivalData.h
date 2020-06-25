//
// Created by lostrong on 18/7/1.
//

#ifndef TOAIN_ARRIVALDATA_H
#define TOAIN_ARRIVALDATA_H
#include<vector>

class ArrivalData{
    public:
    vector<std::pair<double, int>> full_list;
    vector<int> arrival_nodes;
    long total_queries_plan;
    long total_updates_plan;
        ArrivalData(vector <std::pair<double, int>> full_list, vector<int> arrival_nodes, long total_queries_plan,
        long total_updates_plan){
            this->arrival_nodes.insert(this->arrival_nodes.end(), arrival_nodes.begin(), arrival_nodes.end());
            this->full_list.insert(this->full_list.end(), full_list.begin(), full_list.end());
            this->total_queries_plan=total_queries_plan;
            this->total_updates_plan=total_updates_plan;
        }

};

#endif //TOAIN_ARRIVALDATA_H
