//
// Created by lostrong on 18/3/10.
//

#ifndef TOAIN_MULTICOREVARIABLES_H
#define TOAIN_MULTICOREVARIABLES_H

#include <sys/time.h>
#include "../graph/CellNode.h"
#include "../../VTree-master/vtree_knn_demo.h"

vector<KNode> merge_queue;

#define VERIFY 0

#define TOAIN_UPDATE 1

#define TOAIN_QUERY 0

int is_simulation = 0;

int overload_flag = 0;

unordered_map<string, string> paras;

MultiTestParameter multiTestPara;

vector<long> response_time_list;

vector<long> update_response_time_list;

vector<long> query_time_list;

vector<long> update_time_list;

//int is_thresholded=0;

//int is_single_aggregate=0;

const int MAX_CORES = 40;

vector<vector<KNode> > verify_results;

vector<G_Tree*> vtrees;

int** vtree_objects;

//string method_name="dijk";

#define QUERY_ID_FOLD 20000000

vector<KNode> **hier_local_knn_arr_multi;

double update_finish_rate=0.0;

double query_finish_rate = 0.0;

const int rand_length=50000000;

int global_random_numbers[rand_length*2];

//int num_threads_update;

int num_threads_query;
struct timeval global_start;

long long number_of_queries=1;
//double total_query_time=0.0; // secs

int* car_id_insert_cnt;

int can_estimate = 0;

int not_record = 0;

int query_gap=1;

//int query_node = 0;

class GlobalThreadVar {
public:

    int stop_run_threads = 0;
    std::mutex ran_global_locker;
    vector<int> ran_threshold;
    double total_query_time=0.0;
    int number_of_queries=0;


    GlobalThreadVar(){
        stop_run_threads=0;
        total_query_time=0.0;
        number_of_queries=0;
    }
};

GlobalThreadVar **globalThreadVar[2];
long number_of_updates=1;
double total_update_response_time=0.0;
double total_update_process_time=0.0;
long number_of_query_processings=1;
double total_query_process_time = 0.0;
double avg_offset=0.0;
std::mutex update_time_mutex;
std::mutex estimate_mutex;

std::mutex vtree_read_mutex;

#endif //TOAIN_MULTICOREVARIABLES_H
