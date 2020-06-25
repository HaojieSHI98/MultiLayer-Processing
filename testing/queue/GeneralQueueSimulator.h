//
// Created by lostrong on 17/3/9.
//

#ifndef SOB_SIMULATORFORTHROUGHPUT_H
#define SOB_SIMULATORFORTHROUGHPUT_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <queue>

using namespace std;

#define QUERY 0
#define UPDATE 1
#define INSERT 1
#define DELETE 2

int process_time[3] = {58, 6, 6};  // 0: query; 1: insert; 2: delete

struct Query{
    int node_id;
    double init_time; // in seconds
    int operator_type;  // 0: QUERY; 1: UPDATE
    int process_time;  // in microsecs
    Query(int id, double t1, int type, int t2) {
        node_id = id;
        init_time = t1;
        operator_type = type;
        process_time = t2;
    }
};

float nextTime(double rateParameter) {
    return -log(1.0f - (double) rand() / (RAND_MAX + 1.0)) / rateParameter;
}

vector<double> poisson(double lambda, int T) {// T: seconds
    vector<double> sequence;
    double t=0;
//    srand(time(NULL));
    int num=0;
    while(t < T) {
        double element = nextTime(lambda);
        if(t+element < T) {
            sequence.push_back(t + element);
            //cout<< t+element<<endl;
        }
        t+=element;
        num++;
    }
    return sequence;
}

// T: seconds
vector<Query> generate_queryList(double lambda, int T, vector<int> & process_time_list, int type) {
    vector<Query> list;
    vector<double> sequence = poisson(lambda, T);
//    cout<<endl<<sequence.size()<<" "<<lambda*T<<endl;
    int m=process_time_list.size();
    int n = sequence.size();
    for(int i=0; i < n; i++) {
        int id = rand()%n;
        Query query(id, sequence[i], type, process_time_list[i%m]);
        list.push_back(query);
    }
    return list;
}

Query* getNext(vector<Query> &queryList, vector<Query> &insertList, vector<Query> &deleteList, int &i, int &j, int &k) {

    int kind=-1;
    int index=-1;
    double g_time=-1;

    if( i < queryList.size() ) {
        kind = 0;
        index= i;
        g_time = queryList[i].init_time;
    }

    if(j < insertList.size() ){
        if(kind == -1 || insertList[j].init_time < g_time)
        {
            kind = 1;
            index= j;
            g_time = insertList[j].init_time;
        }
    }

    if(k < deleteList.size()) {
        if(kind == -1 || deleteList[k].init_time < g_time)
        {
            kind = 2;
            index= k;
            g_time = deleteList[k].init_time;
        }
    }

    if(kind == -1) return NULL;
    else
    {
        if(kind == 0) {
            i++;
            return &queryList[index];
        }else if(kind == 1){
            j++;
            return &insertList[index];
        }else {
            k++;
            return &deleteList[index];
        }
    }
}

pair<double, double> simulator_FIFO(vector<double>& simulated_query_costs, vector<double>& simulated_update_costs, vector<Query> &queryList, vector<Query> &insertList, vector<Query> &deleteList, int T){
    assert(T < 10000);
    simulated_query_costs.clear();
    simulated_update_costs.clear();
    int total_time = T * 1000000;
    int i=0, j=0, k=0;
    int count = 0; // the number of queries processed.

    double avg_response_time = 0.0;
    Query* current_query;

    double c_time=0;  // current time for processing queries

    while(c_time <= total_time) {
        current_query = getNext(queryList, insertList, deleteList, i, j, k);
        if (current_query == NULL) {
            break;
        }
        c_time = max(c_time, current_query->init_time * 1000000.0);
        if(c_time + current_query->process_time <= total_time) {
            if(current_query->operator_type == QUERY) {
                count++;
                simulated_query_costs.push_back(current_query->process_time);
                avg_response_time +=  c_time -  current_query->init_time * 1000000 + current_query->process_time;
            }
            else{
                simulated_update_costs.push_back(current_query->process_time);
            }
            c_time +=  current_query->process_time;
        } else {
            break;
        }
    }
    //cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
    if(simulated_update_costs.size()<simulated_update_costs.size()*0.95)
        return pair<double,double>(count * 1.0 / T, -1);
    else
        return pair<double,double>(count * 1.0 / T, avg_response_time / count);
}

pair<double, double> simulator_FIFO_divided_by_updates(vector<double>& simulated_query_costs, vector<double>& simulated_update_costs,
                                                       vector<Query> &queryList, vector<Query> &insertList, vector<Query> &deleteList, int T,
                                                       vector<vector<int> > &query_cost, vector<int> &update_cost){
    int microsecs_per_sec=1000000;
    assert(T < 10000);
    simulated_query_costs.clear();
    simulated_update_costs.clear();
    int total_time = T * microsecs_per_sec;
    int i=0, j=0, k=0;
    int count = 0; // the number of queries processed.

    double avg_response_time = 0.0;
    Query* current_query;

    double c_time=0;  // current time for processing queries
    int query_slot;
    int query_processed_in_one_update_slot = 0;
    while(c_time <= total_time) {
        current_query = getNext(queryList, insertList, deleteList, i, j, k);
        if (current_query == NULL) {
            break;
        }
        query_slot = j % query_cost.size();
        if(current_query->operator_type == QUERY)
             current_query->process_time = query_cost[query_slot][query_processed_in_one_update_slot % query_cost[query_slot].size()];
        else
            current_query->process_time =update_cost[query_slot];

            c_time = max(c_time, current_query->init_time * microsecs_per_sec);
//        if(c_time + current_query->process_time <= total_time) {
        if(c_time<=total_time){
            if(current_query->operator_type == QUERY) {
                simulated_query_costs.push_back(current_query->process_time);
                count++;
                query_processed_in_one_update_slot++;
                avg_response_time +=  c_time -  current_query->init_time * microsecs_per_sec + current_query->process_time;

            }else
            {
                simulated_update_costs.push_back(current_query->process_time);
                query_processed_in_one_update_slot = 0;

            }
            c_time +=  current_query->process_time;

        } else {
            break;
        }
    }
//    cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
//    if(simulated_update_costs.size()<insertList.size()*0.95)// the condition used for the current version of paper
    if((simulated_update_costs.size()+simulated_query_costs.size())<(insertList.size()+queryList.size())*0.99)
        // then not satisfy the conditions, because the simulated sizes are less than 99% of the task size.
        return pair<double,double>(count * 1.0 / T, -1);
        // cannot handle all the queries and updates
    else
        return pair<double,double>(count * 1.0 / T, avg_response_time / count);
}

pair<double, double> simulator_QueryFirst(vector<double>& simulated_query_costs, vector<double>& simulated_update_costs,vector<Query> &queryList, int T, double period_time, vector<int> &batch_update_costs){
    assert(T < 10000);
    simulated_query_costs.clear();
    simulated_update_costs.clear();
    int total_time; // terminal time for each period
    int count = 0; // the number of queries processed.

    double avg_response_time = 0.0;
    Query* current_query;

    double c_time=0;  // current time for processing queries
    int microsecs_per_sec=1000000;
    T = T*microsecs_per_sec; // transform T to microseconds
    //int num_updates = T/period_time;
    int current_update_index=0;

    assert(batch_update_costs.size()>0);

    bool finished = false;
    int update_time_rest;

    while(current_update_index * period_time < T){

        c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];

        total_time = min((int)period_time * (current_update_index+1), T);

        update_time_rest = batch_update_costs[current_update_index % batch_update_costs.size()];
        simulated_update_costs.push_back(update_time_rest);
        while(c_time + update_time_rest <= total_time) {
            if( count < queryList.size()) {
                current_query = &queryList[count];
            }else {
                finished = true;
                break;
            }

            if (current_query->init_time*1000000.0 > c_time){
                update_time_rest -= (current_query->init_time*1000000.0 - c_time);
                c_time = current_query->init_time*1000000.0;
            }

            //c_time = max(c_time, current_query->init_time*1000000.0);
            if(c_time + current_query->process_time +update_time_rest <= total_time) {
                count++;
                simulated_query_costs.push_back(current_query->process_time);
                avg_response_time +=  c_time -  current_query->init_time * 1000000 + current_query->process_time;
                c_time += current_query->process_time;
            } else {
                current_query->process_time = current_query->process_time - (total_time - c_time - update_time_rest);
                break;
            }
        }
        if(finished) break;
        current_update_index++;
    }

    //cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
    return pair<double,double>(count * 1.0 / T*microsecs_per_sec, avg_response_time / count);
}

pair<double, double> simulator_QueryFirst_divided_by_updates(vector<double> & simulated_query_costs, vector<double> & simulated_update_costs,
                                                             vector<Query> &queryList, int T, double period_time, vector<int> &batch_update_costs, vector<vector<int> > &query_cost){
    simulated_query_costs.clear();
    simulated_update_costs.clear();
    assert(T < 10000);
    int total_time; // terminal time for each period
    int count = 0; // the number of queries processed.

    double avg_response_time = 0.0;
    Query* current_query;

    double c_time=0;  // current time for processing queries
    int microsecs_per_sec=1000000;
    T = T*microsecs_per_sec; // transform T to microseconds
    int num_updates = T/period_time;
    int current_update_index=0;

    assert(batch_update_costs.size()>0);

    bool finished = false;

    int query_processed_in_one_update_slot = 0;

    int update_time_rest;

    while(current_update_index * period_time < T){

//        cout<<current_update_index<<endl;

        c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];

        total_time = min((int)period_time * (current_update_index+1), T);

        update_time_rest = batch_update_costs[current_update_index % batch_update_costs.size()];

        simulated_update_costs.push_back(update_time_rest);

        query_processed_in_one_update_slot = 0;

        while(c_time + update_time_rest <= total_time) {
            if( count < queryList.size()) {
                current_query = &queryList[count];
            }else {
                finished = true;
                break;
            }

            current_query->process_time = query_cost[current_update_index % query_cost.size()][query_processed_in_one_update_slot % query_cost[current_update_index % query_cost.size()].size()];

            if (current_query->init_time*1000000.0 > c_time){
                update_time_rest -= (current_query->init_time*1000000.0 - c_time);
                c_time = current_query->init_time*1000000.0;
            }

            //c_time = max(c_time, current_query->init_time*1000000.0);
            if(c_time + current_query->process_time + update_time_rest <= total_time) {
                count++;
                query_processed_in_one_update_slot++;
                avg_response_time +=  c_time -  current_query->init_time * 1000000 + current_query->process_time;
                simulated_query_costs.push_back(current_query->process_time);
                c_time += current_query->process_time;
            } else {
                current_query->process_time = current_query->process_time - (total_time - c_time - update_time_rest);
                break;
            }
        }
        if(finished) break;
        current_update_index++;
    }


    //cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
    return pair<double,double>(count * 1.0 / T*microsecs_per_sec, avg_response_time / count);
}

/**
 * Function: simulate the case of "uber's model", that is, the objects are updated batch by batch, and the query can
 * preempt the updates. Every time $period_time, the objects are updated in a batch (the total cost is a value in $batch_update_costs).
 * The queries are executed one by one, whose costs are recorded in query_costs. Note that, the lengths of $query_costs and $batch_update_costs
 * should be enough for the simulation for a period of $simulation_time.
 *
 * Procedure: within each period, the simulation performs like this: 1) pick a cost from $batch_update_costs, subtract from the time cost for update. 2)
 * in the remaining time, do simulations for the queries. 3) calculate the throughput (per sec).
 *
 * @param threshold_time   the threshold for the expected response time (in microsecs)
 * @param period_time    the periodical batch updates of objects (in microsecs)
 * @param query_costs    the list of query costs (in microsecs), where query_costs[i] contains the costs after the i-th update
 * @param batch_update_costs   the list of batch update costs (in microsecs)
 * @param simulation_time    the duration of doing simulation (in secs)
 * @return throughput
 */
int test_for_queryfirst_divided_by_updates(vector<double> & simulated_query_costs, vector<double>& simulated_update_costs,
                                           int simulation_time, double threshold_time, double period_time, vector<vector<int> >& query_costs, vector<int>& batch_update_costs){

    double l = 0, r = 100000; // be careful of 'r'
    int throughput = 0;
    int lambda;
    while(l <= r){
        lambda = (l+r)/2;
//        cout<<"lambda: "<<lambda<<endl;
        vector<Query> queryList = generate_queryList(lambda, simulation_time, batch_update_costs, QUERY);  // ignore "batch_update_cost"
//        cout<<"the size of querylist: "<<queryList.size()<<endl;
        pair<double,double> result = simulator_QueryFirst_divided_by_updates(simulated_query_costs, simulated_update_costs,
                                                                             queryList, simulation_time, period_time, batch_update_costs, query_costs);

        if(result.second <= threshold_time){
            if(throughput < result.first) {
//                cout<<"here:" << result.first<<endl;
                throughput=result.first;
            }
            l = lambda + 10;
//            l=l*1.1;
        }else{
            r = lambda - 10;
//            r=r/1.1;
        }
        cout<<lambda<<endl;
    }
    return throughput;

}

/**
 * Function: simulate the case of "uber's model", that is, the objects are updated batch by batch, and the query can
 * preempt the updates. Every time $period_time, the objects are updated in a batch (the total cost is a value in $batch_update_costs).
 * The queries are executed one by one, whose costs are recorded in query_costs. Note that, the lengths of $query_costs and $batch_update_costs
 * should be enough for the simulation for a period of $simulation_time.
 *
 * Procedure: within each period, the simulation performs like this: 1) pick a cost from $batch_update_costs, subtract from the time cost for update. 2)
 * in the remaining time, do simulations for the queries. 3) calculate the throughput (per sec).
 *
 * @param threshold_time   the threshold for the expected response time (in microsecs)
 * @param period_time    the periodical batch updates of objects (in microsecs)
 * @param query_costs    the list of query costs (in microsecs)
 * @param batch_update_costs   the list of batch update costs (in microsecs)
 * @param simulation_time    the duration of doing simulation (in secs)
 * @return throughput
 */
int test_for_queryfirst(vector<double>& simulated_query_costs, vector<double>& simulated_update_costs, int simulation_time, double threshold_time, double period_time, vector<int>& query_costs, vector<int>& batch_update_costs){
    int l = 0, r = 10000000; // be careful of 'r'
    int throughput = 0;
    int lambda;
    while(l <= r){
        lambda = (l+r)/2;
        vector<Query> queryList = generate_queryList(lambda, simulation_time, query_costs, QUERY);
        pair<double,double> result = simulator_QueryFirst(simulated_query_costs, simulated_update_costs, queryList, simulation_time, period_time, batch_update_costs);
        if(result.second <= threshold_time){
            if(throughput < result.first) throughput=result.first;
            l = lambda + 1;
        }else{
            r = lambda - 1;
        }
    }
    return throughput;
}


/**
 * Function: simulate the case of "Pokemon's model", that is, the objects are updated individually, and the priority of the queries
 * and updates are the same.
 * @param simulation_time the duration of doing simulation (in secs)
 * @param threshold_time  the threshold for the expected response time (in microsecs)
 * @param obj_update_ratio  the frequency of the object updates per second, corresponding to the lambda in Possion process.
 * @param query_costs   the list of query costs (in microsecs), where the costs in query_costs[i] are divided by i-th update
 * @param update_costs   the list of update costs
 * @return throughput
 */
int test_for_fifo_divided_by_updates(vector<double>& simulated_query_costs, vector<double>& simulated_update_costs, int simulation_time, double threshold_time, double obj_update_ratio, vector<vector<int> >& query_costs, vector<int>& update_costs){

    int l = 0, r = 200000; // be careful of 'r'
    int throughput = 0;
    int lambda;
    while(l <= r){
        lambda = (l+r)/2;
        vector<Query> queryList = generate_queryList(lambda, simulation_time, update_costs, QUERY); //ignore "update_costs"

//        cout<<"query list size: "<<queryList.size()<<endl;
//        cout<<"lambda: "<<lambda<<endl;
        vector<Query> insertList = generate_queryList(obj_update_ratio, simulation_time, update_costs, UPDATE);
//        cout<<"insert list size: "<<insertList.size()<<endl;
        vector<Query> deleteList; //empty
        pair<double,double> result = simulator_FIFO_divided_by_updates(simulated_query_costs, simulated_update_costs,
                                                                       queryList, insertList, deleteList, simulation_time, query_costs, update_costs);

        if(result.second>0 && result.second <= threshold_time){
            if(throughput < result.first) throughput=result.first;
            l = lambda + 10;
        }else{
            r = lambda - 10;
        }
    }
    return throughput;
}

/**
 * Function: simulate the case of "Pokemon's model", that is, the objects are updated individually, and the priority of the queries
 * and updates are the same.
 * @param simulation_time the duration of doing simulation (in secs)
 * @param threshold_time  the threshold for the expected response time (in microsecs)
 * @param obj_update_ratio  the frequency of the object updates per second, corresponding to the lambda in Possion process.
 * @param query_costs   the list of query costs
 * @param update_costs   the list of update costs
 * @return throughput
 */
int test_for_fifo(vector<double>& simulated_query_costs, vector<double>& simulated_update_costs, int simulation_time, double threshold_time, double obj_update_ratio, vector<int>& query_costs, vector<int>& update_costs){
    int l = 0, r = 1000000; // be careful of 'r'
    int throughput = 0;
    int lambda;
    while(l <= r){
        lambda = (l+r)/2;
        vector<Query> queryList = generate_queryList(lambda, simulation_time, query_costs, QUERY);
        vector<Query> insertList = generate_queryList(obj_update_ratio, simulation_time, update_costs, UPDATE);
        vector<Query> deleteList; //empty
        pair<double,double> result = simulator_FIFO(simulated_query_costs, simulated_update_costs, queryList, insertList, deleteList, simulation_time);

        if(result.second>0 && result.second <= threshold_time){
            if(throughput < result.first) throughput=result.first;
            l = lambda + 1;
        }else{
            r = lambda - 1;
        }
    }
    return throughput;
}
int main_test() {
    //int T=10;
    // double lambda=15486;
    // vector<Query> queryList = generate_queryList(lambda, T, 0);
    // vector<Query> insertList = generate_queryList(lambda, T, 1);
    // vector<Query> deleteList = generate_queryList(lambda, T, 2);
    //double throughput = simulator(queryList, insertList, deleteList, T, 0);
    //cout<< "throughput: "<< throughput << "  "<< ( 1000000 - lambda * 6 )/ 58.0 << endl;
    return 0;
}
#endif //SOB_SIMULATORFORTHROUGHPUT_H
