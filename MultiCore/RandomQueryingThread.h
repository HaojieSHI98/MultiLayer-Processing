//
// Created by lostrong on 18/3/10.
//

#ifndef TOAIN_QUERYINGTHREAD_H
#define TOAIN_QUERYINGTHREAD_H

//linux, g++ -std=c++14 -o t *.cpp -pthread
#include <queue>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>
#include "../MultiCore/DijkstraKNN.h"
#include "../util/MemoryAllocation.h"
#include "../graph/CellNode.h"
#include "../queue/GeneralQueueSimulator.h"
#include "../MultiCore/Simulation.h"
#include "../MultiCore/ComputeKStar.h"
#include "../MultiCore/MultiCoreVariables.h"
#include "../MultiCore/FunctionTemplate.h"
#include <sys/time.h>
#include "util.h"
#include <condition_variable>
#include<numeric>

//std::mutex global_locker;


class RandomAggregateThread {
private:
    int copy_id;
    std::thread _thread;
    bool is_free;
    queue<pair<int, vector<KNode> > > partial_result_queue;//first query_id; second top-k(k_star) list
    queue<long> issue_time_queue;
    int num_threads_update;
    vector<KNode> knn_result[QUERY_ID_FOLD];
    vector<int> merge_cnt;
    int thread_stop;
    int k;
    std::mutex thread_mutex;
    timeval end;
    std::condition_variable cv;
    bool is_wait=false;
    int pool_id;



public:

    //构造
    RandomAggregateThread(int pool_id_val,int copy_id_val, int k_val, int num_threads_update_val) : is_free(true), thread_stop(0) {
        copy_id = copy_id_val;
        k = k_val;
        pool_id = pool_id_val;
        num_threads_update = num_threads_update_val;
        for (int i = 0; i < QUERY_ID_FOLD; i++) merge_cnt.push_back(0);
        _thread = std::thread(&RandomAggregateThread::run, this);
        // working_thread.detach(); //放到后台， join是等待线程结束
    }

    void join() {
        if(_thread.joinable())
            cout<<"join!"<<endl;
            _thread.join();
            cout<<"ok!"<<endl;
    }

    void set_stop() {
        thread_stop = 1;
    }

    vector<KNode> merge_k(vector<KNode> &S1, vector<KNode> &S2, int k) {
//        cout<<"merge_k"<<endl;
        vector<KNode> tmp;
        int i = 0;
        int j = 0;
        int cnt = 0;
        while (cnt < k) {
            if (i < S1.size() && j < S2.size()) {
                if (S1[i].dis < S2[j].dis) {
                    tmp.push_back(S1[i]);
                    i++;
                } else {
                    tmp.push_back(S2[j]);
                    j++;
                }
                cnt++;
            } else if (i < S1.size()) {
                tmp.push_back(S1[i]);
                i++;
                cnt++;
            } else if (j < S2.size()) {
                tmp.push_back(S2[j]);
                j++;
                cnt++;
            } else {
                break;
            }

        }
        return tmp;
    }

    void add_task(long issue_time, pair<int, vector<KNode> > &partial_result) {
        thread_mutex.lock();
        issue_time_queue.push(issue_time);
        partial_result_queue.push(partial_result);

//        if(is_wait)
            cv.notify_all();
        thread_mutex.unlock();
//        cout<<"aggregate task size: "<<partial_result_queue.size()<<endl;
    }

    void notify(){
        cv.notify_all();
    }

    void run() {
        if(DISPLAY)
        {
            cout << "enter aggregate run!" << endl;
        }
        while (true) {
            if(overload_flag) break;
            if(thread_stop){
                break;
            }
//            if(!thread_stop)
//             {
                 if(!is_wait) {
                     if (partial_result_queue.empty()) {
                         std::unique_lock<std::mutex> aggregate_lock(thread_mutex);
//                         cv.wait(aggregate_lock, []{return true;});
                         is_wait = true;
                         cv.wait(aggregate_lock, [this] { return !(partial_result_queue.empty()) || thread_stop; });
                         is_wait = false;
                     }
                 }


//            }

            if (!partial_result_queue.empty()) {

                thread_mutex.lock();
                // DO NOT use reference here!!
                pair<int, vector<KNode> > partial_result = partial_result_queue.front();
                long issue_time = issue_time_queue.front();
                partial_result_queue.pop();
                issue_time_queue.pop();

                if(can_estimate)
                    gettimeofday(&end, NULL);
                else{
                    estimate_mutex.lock();
                    gettimeofday(&end, NULL);
                    estimate_mutex.unlock();
                }
                long current_time1 =
                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                int adjust_id = partial_result.first % QUERY_ID_FOLD;
                merge_cnt[adjust_id]++;
                if (merge_cnt[adjust_id] == 1) {
                    knn_result[adjust_id]=partial_result.second;
//                    cout<<"yyyyyyyyyyyyyy: "<<adjust_id<<" "<<knn_result[adjust_id].size()<<endl;
                } else {
//                    cout<<"xxxxxxxxxxxxx!"<<endl;
                    knn_result[adjust_id] = merge_k(knn_result[adjust_id], partial_result.second, k);
                }
                thread_mutex.unlock();

                if (merge_cnt[adjust_id] == num_threads_update) {

                    if(can_estimate)
                        gettimeofday(&end, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&end, NULL);
                        estimate_mutex.unlock();
                    }
                    thread_mutex.lock();
                    long current_time =
                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
                    long response_time = current_time - issue_time;
//                    thread_mutex.lock();
//                    response_time_list.push_back(response_time);
                    thread_mutex.unlock();
//                    if(response_time>0.1*MICROSEC_PER_SEC){
//                        cout<<"pool "<<pool_id<<"response time too large:"<<response_time<<endl;
////                        cout<<"overloaded!"<<endl;
////                        overload_flag = 1;
////                        break;
//
//                    }
//                    cout<<"current_time: "<<current_time<<endl;
//                    cout<<"issue_time: "<<issue_time<<endl;
                    if(No_Record_Flag==0)
                    {
                        globalThreadVar[pool_id][copy_id]->ran_global_locker.lock();
                        globalThreadVar[pool_id][copy_id]->total_query_time += response_time * 1.0 / MICROSEC_PER_SEC;
                        globalThreadVar[pool_id][copy_id]->number_of_queries++;
                        globalThreadVar[pool_id][copy_id]->ran_global_locker.unlock();
                    }


//                    observer.tq_mutex[pool_id].lock();
//                    if(observer.tq_ex[pool_id].size()>EXP_SIZE)
//                    {
//                        observer.tq_ex[pool_id].erase(observer.tq_ex[pool_id].begin(),observer.tq_ex[pool_id].begin()+1);
//                    }
//                    observer.tq_ex[pool_id].push_back(response_time);
//                    observer.tq_mutex[pool_id].unlock();
                    if(X_STAR_MODE)
                    {
                        observer.ta_mutex[pool_id].lock();
                        if(observer.ta[pool_id].size()>EXP_SIZE)
                        {
                            observer.ta[pool_id].erase(observer.ta[pool_id].begin(),observer.ta[pool_id].begin()+1);
                        }
                        observer.ta[pool_id].push_back(current_time-current_time1);
                        observer.ta_mutex[pool_id].unlock();
                    }

                    // cout top-k
                    thread_mutex.lock();
                    merge_cnt[adjust_id] = 0;
                    thread_mutex.unlock();
                    if(VERIFY) {
                        cout<<"start verifying"<<endl;
                        int res = verify_two_results_simple(knn_result[adjust_id], verify_results[adjust_id], k);
                        if(res){
                            cout<<"correct!"<<endl;
                        }
                    }

                    if (multiTestPara.is_thresholded) {
                        globalThreadVar[pool_id][copy_id]->ran_global_locker.lock();
                        globalThreadVar[pool_id][copy_id]->ran_threshold[adjust_id] = INT_MAX;
                        globalThreadVar[pool_id][copy_id]->ran_global_locker.unlock();
                    }
                }

            }

//            else if (thread_stop) {
            thread_mutex.lock();
            bool stop_cond = partial_result_queue.empty() && thread_stop;
            thread_mutex.unlock();
            if (stop_cond) {
                if(!multiTestPara.is_single_aggregate) {
                    if (globalThreadVar[pool_id][copy_id]->stop_run_threads == num_threads_update) {
                        cout << "stopped aggregate threads running" << endl;
                        break;
                    }
                }
                else{
                    if (globalThreadVar[pool_id][copy_id]->stop_run_threads == num_threads_update*multiTestPara.num_threads_query) {
                        cout << "stopped aggregate threads running" << endl;
                        break;
                    }

                }
            }
        }
//        cout << endl << endl << "out aggregate" << endl;
    }

};

class RandomThread {
private:
    int copy_id;
    std::thread working_thread;
    int thread_id;
    bool is_free;
    bool is_wait;
    queue<pair<long, pair<int, int> > > task_queue;
    int task_size = 0;
    vector<pair<long, pair<int, int> > > _task_cache_array;
    int _task_index=0;
    int thread_stop;
    int k;
    mem_struct mems;
    std::mutex thread_mutex;
    //Algorithm data structure
    int *dijkstra_object_map;
    int query_id = 0;
    RandomAggregateThread *aggregate_thread;
    std::condition_variable cv;
    timeval end;
    vector<KNode>* hier_local_knn_arr_single;
    int last_query_cost;
    int last_insert_cost;
    int last_delete_cost;
    long last_query_response_time;
    int num_queries_in_queue;
    int num_inserts_in_queue;
    int num_deletes_in_queue;
    int pool_id;


public:

    //构造
    RandomThread(int pool_id_val, int copy_id_val, int id, int k_val,int query_cost,int insert_cost,int delete_cost, RandomAggregateThread *aggregate_thread_val) : thread_id(id), is_free(true),
                                                                                                    thread_stop(0) {
        copy_id = copy_id_val;
        pool_id = pool_id_val;
        k = k_val;
        is_wait=false;
        allocate_mem(mems, test_n + test_n + 1);
        dijkstra_object_map = new int[test_n + 1];
        memset(dijkstra_object_map, 0, sizeof(int) * (test_n + 1));
        aggregate_thread = aggregate_thread_val;

        if(multiTestPara.method_name.compare("toain")==0){
            hier_local_knn_arr_single= new vector<KNode>[test_n+1];
        }
        num_queries_in_queue=0;
        num_inserts_in_queue=0;
        num_deletes_in_queue=0;
        last_query_cost=query_cost;
        last_insert_cost=insert_cost;
        last_delete_cost=delete_cost;
        last_query_response_time=0;
//        cout<<"query_lost:"<<last_query_cost<<"insert_lost:"<<last_insert_cost<<"delete_lost:"<<last_delete_cost<<endl;

        working_thread = std::thread(&RandomThread::run, this);
        // working_thread.detach(); //放到后台， join是等待线程结束
    }

    ~RandomThread(){
        delete_mems(mems);
        if(multiTestPara.method_name.compare("toain")==0) {
            delete[] hier_local_knn_arr_single;
        }
    }

    int get_num_queries_in_queue(){
        return num_queries_in_queue;
    }

    int get_num_inserts_in_queue(){
        return num_inserts_in_queue;
    }

    int get_num_deletes_in_queue(){
        return num_deletes_in_queue;
    }

    void join() {
        if(working_thread.joinable())
            working_thread.join();
    }

    void set_stop() {
        thread_stop = 1;
    }

    long get_est_cost(){
//        return last_query_response_time;
        return last_insert_cost * num_inserts_in_queue + last_delete_cost * num_deletes_in_queue + last_query_cost * num_queries_in_queue;
    }
    void add_task(pair<long, pair<int, int> >& taskids) {
//        if (is_free) {
        // std::lock_guard<std::mutex>(_locker2);
        thread_mutex.lock();
        task_queue.push(taskids);
        task_size++;
//        cout<<"query task size: "<<num_queries_in_queue<<" update task size: "<<num_inserts_in_queue+num_deletes_in_queue
//            << " thread id: "<<thread_id<<" copyid: "<<copy_id<< endl;
        //taskids.second.second is the type of the task
        if(taskids.second.second == QUERY){
            num_queries_in_queue++;
        }
        if(taskids.second.second == INSERT){
            num_inserts_in_queue++;
        }
        if(taskids.second.second == DELETE){
            num_deletes_in_queue++;
        }
//            is_free = false;
        thread_mutex.unlock();
//        }
//        if(is_wait)
            cv.notify_all();
    }

    void add_task_cache(pair<long, pair<int, int> > taskids) {
        _task_cache_array.push_back(taskids);
        cv.notify_all();

    }

    bool is_task_empty(){
        return _task_index>=_task_cache_array.size();
    }
    int get_task_size(){
        return task_size;
    }

    void notify(){
        cv.notify_all();
    }
    void run() {
        if(DISPLAY)cout << "enter thread run!" << endl;
        timeval init;
        timeval end;
        while (true) {
            // cout<<"start ";
            if(overload_flag) {
                break;
            }
            if(is_simulation && thread_stop){
                break;
            }
//            if(!thread_stop) {
                if(!is_wait) {
                    if (task_queue.empty()) {


                        std::unique_lock<std::mutex> process_lock(thread_mutex);
                        is_wait = true;
                        cv.wait(process_lock, [this] { return !(task_queue.empty()) || thread_stop; });
//                    cv.wait(process_lock, [this] { return true; });
                        is_wait = false;
                    }
                }


//            }

//            if(is_task_empty()){
//                if(!thread_stop) {
//
//                    std::unique_lock<std::mutex> process_lock(thread_mutex);
//                    cv.wait(process_lock, [this] { return !is_task_empty()||thread_stop; });
////                    cv.wait(process_lock, [this] { return true; });
//                }
//
//
//            }
            if (!task_queue.empty()) {

//            if(!is_task_empty()){

                thread_mutex.lock();
                pair<int, int> _taskid = task_queue.front().second;// extract (node, type) pair
                long _task_time = task_queue.front().first;
                task_queue.pop();
                task_size--;

                thread_mutex.unlock();
//                pair<int, int> _taskid = _task_cache_array[_task_index].second;// extract (node, type) pair
//                long _task_time = _task_cache_array[_task_index].first;
//                _task_index++;
                // cout<<" task "<<_taskid.second;

                if (_taskid.second == QUERY) {
//                    cout<<"start querying..."<<endl;
//                    cout<<"query"<<endl;
                    long start = clock();
                    vector<KNode> kNNs;

                    if (multiTestPara.is_thresholded)

                        kNNs = HandleQueryThreshold(copy_id, thread_id, multiTestPara.method_name, k, _taskid.first,
                                                    mems.dist, mems.visited,
                                                    mems.q,
                                                    dijkstra_object_map, globalThreadVar[pool_id][copy_id]->ran_threshold,
                                                    query_id);
                    else{

                        if(can_estimate) {
                            // gettimeofday is not thread safe

                            gettimeofday(&init, NULL);
                        }
                        else{
                            estimate_mutex.lock();
                            gettimeofday(&init, NULL);
                            estimate_mutex.unlock();

                        }

                        kNNs = HandleQuery(copy_id, thread_id, multiTestPara.method_name, k, _taskid.first, mems.dist,
                                           mems.visited, mems.q,
                                           dijkstra_object_map, globalThreadVar[pool_id][copy_id]->ran_threshold, query_id,
                                           hier_local_knn_arr_single);
                        if(can_estimate) {
                            gettimeofday(&end, NULL);
                        }
                        else{
                            estimate_mutex.lock();
                            gettimeofday(&end, NULL);
                            estimate_mutex.unlock();

                        }
                            long processing_time =
                                    (end.tv_sec - init.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - init.tv_usec;
                        long current_time =
                                (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                        last_query_response_time = current_time-_task_time;

                        update_time_mutex.lock();
                        total_query_process_time += processing_time/1000000.0;
//                        query_time_list.push_back(processing_time);
                        number_of_query_processings++;
                        update_time_mutex.unlock();
                        last_query_cost = processing_time;

//                        update_time_mutex.lock();
                        if(X_STAR_MODE)
                        {
                            observer.tq_mutex[pool_id].lock();
                            if(observer.tq_ex[pool_id].size()>EXP_SIZE)
                            {
                                observer.tq_ex[pool_id].erase(observer.tq_ex[pool_id].begin(),observer.tq_ex[pool_id].begin()+1);
                            }
                            observer.tq_ex[pool_id].push_back(processing_time);
                            observer.tq_mutex[pool_id].unlock();
                        }

//                        update_time_mutex.unlock();
                    }
                    if (multiTestPara.is_thresholded) {
                        if (kNNs.size() == k) {
                            globalThreadVar[pool_id][copy_id]->ran_global_locker.lock();

                            if (kNNs[k - 1].dis < globalThreadVar[pool_id][copy_id]->ran_threshold[query_id])
                                globalThreadVar[pool_id][copy_id]->ran_threshold[query_id] = kNNs[k - 1].dis;
                            globalThreadVar[pool_id][copy_id]->ran_global_locker.unlock();
                        }
                    }
//                    if(multiTestPara.num_threads_update>1) {
                        pair<int, vector<KNode> > partial_res = make_pair(query_id, kNNs);
                        query_id = (query_id + 1) % QUERY_ID_FOLD;
                        aggregate_thread->add_task(_task_time, partial_res);
//                    }
//                    else{
//                        gettimeofday(&end, NULL);
//                        long current_time =
//                                (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
//                        long response_time = current_time - _task_time;
//                        if(response_time > 100*MICROSEC_PER_SEC){
//                            cout<<"overloaded!"<<endl;
//                            overload_flag=1;
//                            break;
//                        }
////
//                        globalThreadVar[copy_id]->total_query_time += response_time * 1.0 / MICROSEC_PER_SEC;
//
//                    }
//                    cout<<"aggregate task added!"<<endl;

                    // put partial querying result to somewhere
                }
//                cout<<"method name: "<<multiTestPara.method_name<<endl;
                if (_taskid.second == INSERT) {
//                    cout<<"insert "<<_taskid.first<<endl;

                    if(can_estimate) {
                        gettimeofday(&init, NULL);
                    }
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&init, NULL);
                        estimate_mutex.unlock();

                    }

                    HandleInsert(copy_id, thread_id, multiTestPara.method_name, k, _taskid.first, dijkstra_object_map,
                                 mems,
                                 hier_local_knn_arr_single);
                    if(can_estimate) {
                        gettimeofday(&end, NULL);
                    }
                    else {
                        estimate_mutex.lock();
                        gettimeofday(&end, NULL);
                        estimate_mutex.unlock();
                    }
                    if(not_record==0) {
                        long processing_time =
                                (end.tv_sec - init.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - init.tv_usec;
                        long current_time =
                                (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec -
                                global_start.tv_usec;

                        if (number_of_updates > multiTestPara.init_objects) {
                            update_time_mutex.lock();
                            total_update_response_time += (current_time - _task_time) / 1000000.0;
                            total_update_process_time += processing_time / 1000000.0;
//                            update_time_list.push_back(processing_time);
                            number_of_updates++;
                            update_time_mutex.unlock();

                            if(X_STAR_MODE)
                            {
                                observer.tu_mutex[pool_id].lock();
                                if(observer.tu_ex[pool_id].size()>EXP_SIZE){
                                    observer.tu_ex[pool_id].erase(observer.tu_ex[pool_id].begin(),observer.tu_ex[pool_id].begin()+1);
                                }
                                observer.tu_ex[pool_id].push_back(processing_time);
                                observer.tu_mutex[pool_id].unlock();
                            }
                        }
//                        last_insert_cost = processing_time;
                    }


                }
                if (_taskid.second == DELETE) {
                    // cout<<" delete ";
//                    cout<<"delete "<<_taskid.first<<endl;//

                    if(can_estimate) {
                        gettimeofday(&init, NULL);
                    }
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&init, NULL);
                        estimate_mutex.unlock();

                    }
                    HandleDelete(copy_id, thread_id, multiTestPara.method_name, k, _taskid.first, dijkstra_object_map,
                                 mems, hier_local_knn_arr_single);
                    if(can_estimate) {
                        gettimeofday(&end, NULL);
                    }
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&end, NULL);
                        estimate_mutex.unlock();

                    }
                    if(not_record==0) {
                        long processing_time =
                                (end.tv_sec - init.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - init.tv_usec;
                        long current_time =
                                (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec -
                                global_start.tv_usec;

                        update_time_mutex.lock();
                        total_update_response_time += (current_time - _task_time) / 1000000.0;
                        total_update_process_time += processing_time / 1000000.0;
//                        update_time_list.push_back(processing_time);
                        number_of_updates++;
                        update_time_mutex.unlock();

                        if(X_STAR_MODE)
                        {
                            observer.tu_mutex[pool_id].lock();
                            if(observer.tu_ex[pool_id].size()>EXP_SIZE){
                                observer.tu_ex[pool_id].erase(observer.tu_ex[pool_id].begin(),observer.tu_ex[pool_id].begin()+1);
                            }
                            observer.tu_ex[pool_id].push_back(processing_time);
                            observer.tu_mutex[pool_id].unlock();
                        }
                        // need lock
//                        last_delete_cost = processing_time;
                    }
                    // cout<<" /delete ";

                }

                thread_mutex.lock();
                if(_taskid.second == QUERY){
                    num_queries_in_queue--;
                }
                if(_taskid.second == INSERT){
                    num_inserts_in_queue--;
                }
                if(_taskid.second == DELETE){
                    num_deletes_in_queue--;
                }
                thread_mutex.unlock();

            }

            thread_mutex.lock();
            bool cond_stop = task_queue.empty() && thread_stop;
            thread_mutex.unlock();
//            thread_mutex.lock();
//            bool cond_stop = is_task_empty() && thread_stop;
//            thread_mutex.unlock();
            if (cond_stop) {
//                cout << "stopped running" << endl;

                if(!multiTestPara.is_single_aggregate) {
                    globalThreadVar[pool_id][copy_id]->ran_global_locker.lock();
                    globalThreadVar[pool_id][copy_id]->stop_run_threads++;
                    globalThreadVar[pool_id][copy_id]->ran_global_locker.unlock();
                }
                else {
                    globalThreadVar[pool_id][0]->ran_global_locker.lock();
                    globalThreadVar[pool_id][0]->stop_run_threads++;
                    globalThreadVar[pool_id][0]->ran_global_locker.unlock();
                }
                if(can_estimate) {

                    gettimeofday(&end, NULL);
                }
                else {
                    estimate_mutex.lock();
                    gettimeofday(&end, NULL);
                    estimate_mutex.unlock();
                }
                if(not_record==0) {
                    long current_time =
                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
//                cout<<"copyid: "<<copy_id<<" threadid: "<<thread_id<<" ending time: "<<current_time<<endl;
//                if(multiTestPara.method_name.compare("vtree")==0) {
//                    for (int i = 0; i < test_n; i++) {
//                        if (dijkstra_object_map[i]) {
//                            vtree_objects[copy_id * multiTestPara.num_threads_update + thread_id][i] = 1;
//                        }
//                    }
//                }
                }
                break;
            }
            // cout<<" end"<<endl;
        }
      // cout<<" end"<<endl;
    }
};

typedef struct {
    int threadpool_id;
    int num_thread_update;
    int num_threads_query;
    int _needjoin;
    int run_time;//how many times have run
    int threshold_number;
    vector< vector<int> >total_object_map;
    vector<std::pair<double, int> > init_list;
    vector<int> init_arrival_nodes;
    vector<int> current_object_node;
    int current_frame;//frame
    int begin_frame;
    int global_start_q_id;
    int current_query_threads;
    int rand_idx_query;
    int rand_idx_update;
    long total_offset;
    int total_queries;
    vector<RandomThread *> _pool;
    RandomAggregateThread **_aggregate_thread;
    RandomAggregateThread * _single_aggregate_thread;
    int num_intask;
    double response_time;
    int query_num;
    int restart_flag;
    double response_time_first;
    int query_num_first;
    int ontask_num;
    double eva_response_time;
    int eva_query_num;
    double eva_response_speed;

}OneThreadPool;

class RandomTwoThreadPool_Control {
private:
    std::thread _main_thread;
    int begin_node;
    int end_node;
    long total_queries_plan;
    long total_updates_plan;
    long total_queries_finished;
    long total_updates_finished;
    int* car_nodes;
    string configstr;
    vector<std::pair<double, int> > full_task_list;
    vector<int> arrival_task_nodes;
    int query_cost;
    int insert_cost;
    int delete_cost;
    OneThreadPool tp[2];
    double alpha;
    int test_n;
    double fail_p;
    int k;
    int num_threads_all;
    int num_threads_each;
    vector<std::pair<double, int> > full_list;
    vector<int> arrival_nodes;
    int configurationId;
    int mode;
    int last_eva_t = 0;
    int query_update_set = 1;
    int big_q_query;
    int big_q_update;
    int big_u_query;
    int big_u_update;

public:
    RandomTwoThreadPool_Control( int begin_node_val, int end_node_val, int num_threads_val, double alpha_val, int k_val, double fail_p_val,
                                int test_n_val, string configstr_val,int query_cost_val,int insert_cost_val,int delete_cost_val,int configurationId_val) {
        cout << "constructing RandomThreadPool..." << endl;
        for(int ti = 0;ti<2;ti++)
        {
            tp[ti].threadpool_id=ti;
            tp[ti].run_time = 0;
            tp[ti]._needjoin = 0;
            tp[ti].init_arrival_nodes.clear();
            tp[ti].init_list.clear();
            tp[ti].threshold_number = 0;
            tp[ti].current_object_node.clear();
            tp[ti].total_object_map.clear();
            tp[ti].current_frame = 0;
            tp[ti].begin_frame = 0;
            tp[ti].global_start_q_id = 0;
            tp[ti].num_intask = 0;
            tp[ti].query_num = 0;
            tp[ti].response_time = 0;
            tp[ti].restart_flag = 0;
        }
        mode = NORMAL_MODE;
        tp[1].threshold_number = 20000;
        tp[0].threshold_number = 20000;
        configurationId = configurationId_val;
        begin_node = begin_node_val;
        end_node = end_node_val;
        configstr = configstr_val;
        test_n = test_n_val;
//        if(DISPLAY)cout<<"num_threads_update:"<<num_threads_update<<" num_threads_query:"<<num_threads_query_val<<endl;
        alpha = alpha_val;
        fail_p = fail_p_val;
        k = k_val;
        total_queries_plan=0;
        total_updates_plan=0;
        total_queries_finished=0;
        total_updates_finished=0;
        query_cost = query_cost_val;
        insert_cost = insert_cost_val;
        delete_cost = delete_cost_val;
        num_threads_all = num_threads_val;
        num_threads_each = int((num_threads_all-3)/2);
        first_init();
    }

    //释放线程池
    ~RandomTwoThreadPool_Control() {
        if(DISPLAY)cout << "hi, RandomThreadPooladPool()" << endl;
//         for(int i = 0;i < _pool.size(); ++i){
//             delete _pool[i];
//         }

    }

    vector<int> generate_arrival_nodes(vector<std::pair<double, int> >& full_list, int begin, int end){

        int rand_idx_update=0;
        int rand_idx_query=0;
        vector<int> object_list;
        vector<int> non_object_list;
        vector<int> arrival_node_list;
        for (int i = begin; i < end; i++) {
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
                int index = rand() % non_object_list.size();
                rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;
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
                int index = rand() % object_list.size();
                rand_idx_update=(rand_idx_update+1)%rand_length;
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
                    rand_idx_query=(rand_idx_query+1)%rand_length;
                    query_node = query_node % (test_n-1)+1;

                } while (vStart[query_node] == 0);
                arrival_node_list.push_back(query_node);
            }

        }
        return arrival_node_list;
    }
    void first_init(){
        full_task_list.clear();
        arrival_task_nodes.clear();
        full_list.clear();
        arrival_nodes.clear();
        vector<std::pair<double, int> > full_list;
        int init_objects = multiTestPara.init_objects;
        for (int i = 0; i < init_objects; i++) {
            full_list.push_back(make_pair(0.0, INSERT));
        }
        string myversion = "v1_"+std::to_string(multiTestPara.init_objects)+"_";
        std::ifstream queryfile;

        string queryfile_name = input_parameters.input_data_dir + "query_"+myversion + multiTestPara.configstr+".txt";
        cout<<"query file name:"<<queryfile_name<<endl;
        queryfile.open(queryfile_name,std::ios_base::in);
        double f1;
        int f2;
        if (!queryfile.is_open()) {
            cout << "can't load queryfile!" << endl;
            vector<std::pair<double, int> > append_list = make_online_query_update_list_str(configstr,multiTestPara.init_objects,multiTestPara.layer);
            cout<<"Generated!"<<endl;
            std::ofstream queryfile_w;
            queryfile_w.open(queryfile_name,std::ios_base::out);
            for (pair<double, int> &item : append_list) {
                full_list.push_back(item);
                queryfile_w << item.first << " " << item.second << endl;
            }
            queryfile_w.close();
            cout << "write to queryfile!" << endl;
        } else {
            while (!queryfile.eof()) {
                queryfile >> f1 >> f2;
//                cout<<f1<<" "<<f2<<endl;
                full_list.push_back(make_pair(f1, f2));
            }
            queryfile.close();
            full_list.pop_back();
            cout << "read from queryfile!" << endl;
        }

//        arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);
//        vector<int> arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);

//        vector<int> arrival_nodes;

        std::ifstream nodefile;
        string nodefile_name = input_parameters.input_data_dir + "node_"+myversion +multiTestPara.configstr+".txt";
        cout<<"node file name:"<<nodefile_name<<endl;
        nodefile.open(nodefile_name, std::ios_base::in);

        if(!nodefile.is_open())
        {
            cout<<"can't load nodefile!"<<endl;
            arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);
            std::ofstream nodesfile_w;
            std::ifstream nodesfile_m;
            nodesfile_w.open(nodefile_name, std::ios_base::out);
             for (int node_i=0;node_i<arrival_nodes.size();node_i++)
             {
                 nodesfile_w<<arrival_nodes[node_i]<<endl;
             }
            nodesfile_w.close();
             cout<<"generated nodefile"<<endl;
            nodesfile_m.open(nodefile_name, std::ios_base::in);
            vector<int> arrival_cm;
            while(!nodesfile_m.eof())
            {
                int f3;
                nodesfile_m>>f3;
//            cout<<f3<<endl;
                arrival_cm.push_back(f3);
            }
            nodesfile_m.close();
            cout<<"size:"<<arrival_nodes.size()<<" - "<<arrival_cm.size()<<endl;
            for(int s_i = 0;s_i<arrival_nodes.size();s_i++)
            {
                if(arrival_nodes[s_i]!=arrival_cm[s_i]) cout<<"not equal! "<<s_i<<" "<<arrival_nodes[s_i]<<" - "<<arrival_cm[s_i]<<endl;
            }
            cout<<"Finish Compare!"<<endl<<endl;
        } else{
            while(!nodefile.eof())
            {
                int f3;
                nodefile>>f3;
//            cout<<f3<<endl;
                arrival_nodes.push_back(f3);
            }
            nodefile.close();
            cout<<"read from nodefile!"<<endl;
//            arrival_nodes.erase(arrival_nodes.end()-5,arrival_nodes.end());
            arrival_nodes.pop_back();
            cout<<"full list size:"<<full_list.size()<<" nodes size:"<<arrival_nodes.size()<<endl;
        }

        full_task_list.assign(full_list.begin()+init_objects,full_list.end());
        arrival_task_nodes.assign(arrival_nodes.begin()+init_objects,arrival_nodes.end());
        vector<std::pair<double, int> > init_list;
        init_list.assign(full_list.begin(),full_list.begin()+init_objects);
        vector<int> init_arrival_node;
        init_arrival_node.assign(arrival_nodes.begin(),arrival_nodes.begin()+init_objects);
        cout << "full_list made..." << endl;
        cout << "full list size: " << full_list.size() << endl;

        int query_n = 0,update_n=0;
        for(int f_i=0;f_i<full_task_list.size();f_i++)
        {
            if(full_task_list[f_i].second==QUERY) query_n++;
            else update_n++;
        }
        cout<<"full task list query: "<<query_n<<" update:"<<update_n<<endl<<endl<<endl;
        for(int ti = 0;ti<2;ti++)
        {
            tp[ti].init_list.assign(init_list.begin(),init_list.end());
            tp[ti].init_arrival_nodes.assign(init_arrival_node.begin(),init_arrival_node.end());
        }
//        tp[0].num_threads_query = num_threads_each;
//        tp[0].num_thread_update = 1;
        big_q_query = num_threads_each;
        big_q_update = 1;
        int num_q = int(sqrt(num_threads_each));
        int num_p = int((num_threads_each)/num_q);
        tp[0].num_threads_query = max(num_q,num_p);
        tp[0].num_thread_update = num_p+num_q-max(num_q,num_p);
//        tp[1].num_thread_update = max(num_q,num_p);
//        tp[1].num_threads_query = num_p+num_q-max(num_q,num_p);
        big_u_query = num_p+num_q-max(num_q,num_p);
        big_u_update = max(num_q,num_p);
        tp[1].num_threads_query = num_threads_each;
        tp[1].num_thread_update = 1;
        if(DISPLAY)cout<<"Pool 0 queries:"<<tp[0].num_threads_query<<" updates:"<<tp[0].num_thread_update<<endl;
        if(DISPLAY)cout<<"Pool 1 queries:"<<tp[1].num_threads_query<<" updates:"<<tp[1].num_thread_update<<endl;
        for(int ti=0;ti<2;ti++){
            globalThreadVar[ti] = new GlobalThreadVar*[tp[ti].num_threads_query];
            int k_star = compute_k_star(k, tp[ti].num_thread_update, alpha, fail_p);
            for(int j =0;j<tp[ti].num_threads_query;j++) {
                globalThreadVar[ti][j] = new GlobalThreadVar();
                globalThreadVar[ti][j]->ran_threshold.clear();
//            cout<<"query:"<<j<<" time:"<<globalThreadVar[j]->total_query_time<<" num:"<<globalThreadVar[j]->number_of_queries<<endl;
                for (int i = 0; i <= QUERY_ID_FOLD; i++) globalThreadVar[ti][j]->ran_threshold.push_back(INT_MAX);
            }
            if(!multiTestPara.is_single_aggregate)
                tp[ti]._aggregate_thread = new RandomAggregateThread* [tp[ti].num_threads_query];
            if(multiTestPara.is_single_aggregate){
                tp[ti]._single_aggregate_thread=new RandomAggregateThread(ti,0, k, tp[ti].num_thread_update);
            }
            if(DISPLAY)cout << "k_star: " << k_star << endl;
            for(int j =0;j<tp[ti].num_threads_query;j++){

                if(!multiTestPara.is_single_aggregate)
                    tp[ti]._aggregate_thread[j] = new RandomAggregateThread(ti,j, k, tp[ti].num_thread_update);
                for(int i=0;i<tp[ti].num_thread_update;i++) {
                    if(!multiTestPara.is_single_aggregate) {
                        RandomThread *t = new RandomThread(ti,j, i, k_star,query_cost,insert_cost,delete_cost, tp[ti]._aggregate_thread[j]);
                        tp[ti]._pool.push_back(t);
                    }
                    else{
                        RandomThread *t = new RandomThread(ti,j, i, k_star,query_cost,insert_cost,delete_cost, tp[ti]._single_aggregate_thread);
                        tp[ti]._pool.push_back(t);
                    }
                }
            }
        }
        mem_struct mems;

        if(VERIFY){
            allocate_mem(mems, test_n+1);
            car_nodes = new int[test_n+1];
            memset(car_nodes, 0, sizeof(int)*(test_n+1));
        }
        if(can_estimate)
            gettimeofday(&global_start, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&global_start, NULL);
            estimate_mutex.unlock();
        }
    }
    void join_all(){
        for(int ti = 0;ti<2;ti++){
            join(ti);
        }
        if(_main_thread.joinable())
            _main_thread.join();
        if(DISPLAY)cout << "finish joining main thread" << endl;
    }
    void join(int ti) {
        if(DISPLAY)cout << "start joining" << endl;

        if(!multiTestPara.is_single_aggregate) {
            for (int j = 0; j < tp[ti].num_threads_query; j++) {
                tp[ti]._aggregate_thread[j]->join();
            }
        }
        if(multiTestPara.is_single_aggregate){
            tp[ti]._single_aggregate_thread->join();
        }
        if(DISPLAY)cout << "finish joining aggregatethreads" << endl;



        int max_updates=0;
        for (int i = 0; i < tp[ti]._pool.size(); i++) {
            int remain_updates = tp[ti]._pool[i]->get_num_inserts_in_queue()+tp[ti]._pool[i]->get_num_deletes_in_queue();
            if(remain_updates > max_updates){
                max_updates =  remain_updates;
            }
            tp[ti]._pool[i]->join();
        }

        update_finish_rate = max_updates*1.0/total_updates_plan;
        for (int i = 0; i < tp[ti].num_threads_query; i++) {
            total_queries_finished += globalThreadVar[tp[ti].threadpool_id][i]->number_of_queries;
        }
        query_finish_rate = total_queries_finished * 1.0 / total_queries_plan;

        if(DISPLAY)cout << "finish joining threads" << endl;

    }

    int isNeedJoin() {
        return tp[0]._needjoin&&tp[1]._needjoin;
    }

    void start() {
        _main_thread = std::thread(&RandomTwoThreadPool_Control::run, this);
    }

    void task_init(int id){
        struct timeval end;
        tp[id].current_query_threads=0;
        tp[id].rand_idx_query=0;
        tp[id].rand_idx_update=rand_length;
        tp[id].total_offset=0;
        tp[id].total_queries=0;
        for(int j =0;j<40;j++) {
            vector<int> tmp;
            for (int i = 0; i <= test_n; i++) {
                tmp.push_back(0);
            }
            tp[id].total_object_map.push_back(tmp);
        }
        tp[id].global_start_q_id = 0;
//        cout<<"Pool"<<id<<" step 1 finished!"<<endl;
        for(int init_i=0;init_i<tp[id].init_list.size();init_i++){
            if(can_estimate)
                gettimeofday(&end, NULL);
            else{
                estimate_mutex.lock();
                gettimeofday(&end, NULL);
                estimate_mutex.unlock();
            }
            long current_time =
                    (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

            pair<double, int> &event = tp[id].init_list[init_i];
            int non_object_node = tp[id].init_arrival_nodes[init_i];

            int start_q_id = tp[id].global_start_q_id;
            tp[id].global_start_q_id=(tp[id].global_start_q_id+1)%tp[id].num_thread_update;
            tp[id].rand_idx_update=(tp[id].rand_idx_update+1)%rand_length+rand_length;
//            if(id==1) cout<<"Step 1"<<endl;
            for(int z = 0; z<tp[id].num_threads_query;z++) {


                pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                pair<long, pair<int, int> > task = std::make_pair(current_time, node_type_pair);


                // assign insert tasks
                int min_task_size = INT_MAX;
//                if(id==1) cout<<"start_q_id"<<start_q_id<<"z:"<<z<<endl;
                int pool_index=z * tp[id].num_thread_update + start_q_id;
                tp[id]._pool[pool_index]->add_task(task);
//                if(id==1) cout<<"Step 2"<<endl;
                tp[id].total_object_map[z][non_object_node] = pool_index;
            }
//            if(id==1) cout<<"Step 3"<<endl;
            tp[id].current_object_node.push_back(non_object_node);
            if(VERIFY){
//                    car_nodes[non_object_node]=1;
                DijkstraKNNInsert(non_object_node, car_nodes);
            }
        }
    }
    void update_param(void){
        double t_uq;
        observer.update_rate = 0;
        observer.query_rate = 0;
        for(int o_i=0;o_i<observer.task_list.size();o_i++){
            if(observer.task_list[o_i]==QUERY) observer.query_rate++;
            else observer.update_rate++;
        }
        double time_val = (observer.time_list[observer.task_list.size()-1]-observer.time_list[0])/MICROSEC_PER_SEC;
        observer.update_rate = observer.update_rate/time_val;
        observer.query_rate = observer.query_rate/time_val;
//        observer.last_update_query_ratio = observer.update_query_ratio;
//        observer.task_list.erase(observer.task_list.begin(),observer.task_list.begin()+1);
        if(observer.update_query_ratio_set.size()>=FILTER_NUM)
        {
            observer.update_query_ratio_set.erase(observer.update_query_ratio_set.begin(),observer.update_query_ratio_set.begin()+1);
        }
        if(observer.query_rate<=1e-5) t_uq=100;
        else t_uq = observer.update_rate/observer.query_rate;
        observer.update_query_ratio_set.push_back(t_uq);
        observer.last_update_query_ratio = observer.update_query_ratio;
        if(observer.update_query_ratio_set.size()<FILTER_NUM) observer.update_query_ratio=t_uq;
        else{
            double sum = 0;
            double max = observer.update_query_ratio_set[0];
            double min = max;
            for(int oj =0;oj<observer.update_query_ratio_set.size();oj++)
            {
                sum += observer.update_query_ratio_set[oj];
                if(min>observer.update_query_ratio_set[oj]) min = observer.update_query_ratio_set[oj];
                if(max<observer.update_query_ratio_set[oj]) max = observer.update_query_ratio_set[oj];
            }
            observer.update_query_ratio = (sum-max-min)/(observer.update_query_ratio_set.size()-2);
        }
//        if(observer.query_rate<=1e-5) observer.update_query_ratio=100;
//        else observer.update_query_ratio = observer.update_rate/observer.query_rate;
        if(X_STAR_MODE)
        {
            for(int id = 0;id<2;id++)
            {
                double ta_sum = std::accumulate(observer.ta[id].begin(),observer.ta[id].end(),0);
                double ta_mean = ta_sum/observer.ta[id].size();
                double ts_sum = std::accumulate(observer.ts[id].begin(),observer.ts[id].end(),0);
                double ts_mean = ts_sum/observer.ts[id].size();
                observer.ratio_x[id]=(ta_mean+ts_mean)/tp[id].num_thread_update;
                double tq_sum = std::accumulate(observer.tq_ex[id].begin(),observer.tq_ex[id].end(),0);
                double tq_mean = tq_sum/observer.tq_ex[id].size();
                observer.tq[id] = tq_mean;
                double accum  = 0.0;
                std::for_each (std::begin(observer.tq_ex[id]), std::end(observer.tq_ex[id]), [&](const double d) {
                    accum  += (d-tq_mean)*(d-tq_mean);
                });
                observer.Vq[id] = sqrt(accum/(observer.tq_ex[id].size()-1));

                double tu_sum = std::accumulate(observer.tu_ex[id].begin(),observer.tu_ex[id].end(),0);
                double tu_mean = tu_sum/observer.tu_ex[id].size();
                observer.tu[id] = tu_mean;
                accum  = 0.0;
                std::for_each (std::begin(observer.tu_ex[id]), std::end(observer.tu_ex[id]), [&](const double d) {
                    accum  += (d-tu_mean)*(d-tu_mean);
                });
                observer.Vu[id] = accum/(observer.tu_ex[id].size()-1);

//                cout<<"pool "<<id<<" ta - "<<ta_mean<<" ts - "<<ts_mean<<" tq - "<<tq_mean<<" Vq - "<<observer.Vq[id]<<" tu - "<<tu_mean<<" Vu - "<<observer.Vu[id]<<" Ratio - "<<observer.ratio_x[id]<<endl;
            }
        }


//        cout<<"update_rate:"<<observer.update_rate<<" query_rate:"<<observer.query_rate<<"time_val:"<<time_val<<endl;
    }
    double functionx(int ti,int x,int core){
        double y = floor(core/x);
        double gama_q = observer.Vq[ti]/(observer.tq[ti]*observer.tq[ti]);
        double gama_u = observer.Vu[ti]/(observer.tu[ti]*observer.tu[ti]);
        double t_w = (observer.query_rate*observer.tq[ti]*observer.tq[ti]*(1+gama_q)/y+
                observer.update_rate*observer.tu[ti]*observer.tu[ti]*(1+gama_u)/x)/
                (2*(MICROSEC_PER_SEC-observer.query_rate*observer.tq[ti]/y-observer.update_rate*observer.tu[ti]/x));
        double t_all = t_w+observer.tq[ti]+observer.ratio_x[ti]*x;
//        cout<<"pool: "<< ti<<" x - "<<x<<" tw - "<<t_w<<"t_all - "<<t_all<<" tx - "<<observer.ratio_x[ti]*x<<" tq - "<<observer.tq[ti]<<endl;
        if(t_w<=0) t_all = -1;
        return t_all;
    }
    void compute_x_star(){
        for(int id =0;id<2;id++)
        {
            double t_min=MICROSEC_PER_SEC;
            double x_star=1;
            for(int x = 1;x<=num_threads_each;x++){
                double t_ = functionx(id,x,num_threads_each);
                if(t_<=0) continue;
                
                if(x == 1) t_min = t_;
                if(t_<t_min){
                    x_star = x;
                    t_min = t_;
                }
            }
            if(x_stars.size()>=STAR_NUM)
            {
                x_stars.erase(x_stars.begin(),x_stars.begin()+2);
            }
            x_stars.push_back(x_star);
            if(x_stars.size()>=STAR_NUM)
            {
                for(int i = 1;i<=num_threads_each;i++)
                {
                    int index_num =0;
                    for(int j =1;j<x_stars.size();j++)
                    {
                        if(x_stars[j]==i)
                        {
                            index_num++;
                        }
                    }
                    if(index_num> STAR_NUM/2) cout<<"x_star :"<<i<<" num_of_times: "<<index_num<<endl;
                }
            }
//            cout<<"pool : "<<id<<" x_star : "<<x_star<<endl;
        }

    }
    void task_run(){
        struct timeval end;
        long offset_time;
        int query_turn_flag = 1;
        int turn_num = 0;
        int t_min =1;
        for (int i=0; i < full_task_list.size(); i++) {
            if(overload_flag) break;
            if(arrival_task_nodes[i]==-1) continue;
            pair<double, int> &event = full_task_list[i];

            long issue_time = floor(event.first * MICROSEC_PER_SEC);

//            cout<<"issue_time:"<<issue_time<<endl;
//            int restart_flag = 0;
//            int restart_pool = -1;
//
////            if(event.first>=multiTestPara.config_simulation_time&&tp[0].run_time==0&&tp[1].run_time==0){
//            if(observer.update_rate>20*observer.query_rate&&tp[0].run_time==0&&tp[1].run_time==0&&event.first>=5){
//                update_query_time();
//                clear_query_time();
////                if(tp[0].response_time*tp[1].query_num>tp[1].response_time*tp[0].query_num)
////                    restart_pool = 0;
////                else restart_pool = 1;
//                tp[restart_pool].restart_flag = 1;
//                cout<<"time now: "<<event.first<<endl;
//                cout<<"restart pool"<<restart_pool<<"!!"<<endl<<endl<<endl;
//                task_reinit(0);
//                restart_flag = 1;
//            }


            if(can_estimate)
                gettimeofday(&end, NULL);
            else{
                estimate_mutex.lock();
                gettimeofday(&end, NULL);
                estimate_mutex.unlock();
            }
            long current_time=(end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

            int restart_flag = 0;

            if(event.first>=t_min) {
                update_param();
//                cout<<"uq_r:"<<observer.update_query_ratio<<"uq_r_last:"<<observer.last_update_query_ratio<<endl;
//                cout<<"mode:"<<mode<<" Flag:"<<No_Record_Flag<<endl;
//                for(int id =0;id<2;id++) {
//                    double q_t =0;
//                    int q_n = 0;
//                    for (int i = 0; i < tp[id].num_threads_query; i++) {
//                        q_t += globalThreadVar[id][i]->total_query_time;
//                        q_n += globalThreadVar[id][i]->number_of_queries;
//                    }
//                    cout<<"pool:"<<id<<" response time:"<<q_t<<" num:"<<q_n<<endl;
//                }
                if(mode == NORMAL_MODE && No_Record_Flag == 0)
                {
                    if(i==0) {
                        mode = EVALUATION_START_MODE;
                        last_eva_t = current_time;
                        start_evaluation();
                        cout<<"Last mode: NORMAL--> Start Evaluation!"<<endl<<endl;
//                        No_Record_Flag = 1;
                        Start_No_Record_Time = current_time;
                    }
                    if(abs(observer.update_query_ratio-observer.last_update_query_ratio)>=MIN_UQ_DIFF)
                    {
                        cout<<"changed now!"<<endl;
//                    cout<<"uq_r:"<<observer.update_query_ratio<<"uq_r_last:"<<observer.last_update_query_ratio<<endl;
                        int r_query_num = -1;
                        int r_update_num = -1;
                        if((observer.update_query_ratio-observer.last_update_query_ratio>=MIN_UQ_DIFF)&&(observer.update_query_ratio>=Update_Query_Threshold)){
                            if(tp[0].num_thread_update==big_q_update&&tp[1].num_thread_update==big_q_update)
                            {
                                mode = RESET_MODE;
                                r_query_num = big_u_query;
                                r_update_num = big_u_update;
                            }else if(tp[0].num_thread_update==big_q_update||tp[1].num_thread_update==big_q_update)
                            {
                                mode = EVALUATION_START_MODE;
                                last_eva_t = current_time;
                                start_evaluation();
                                cout<<"Last mode: NORMAL--> Start Evaluation!"<<endl<<endl;
                            }
                        }else if((observer.last_update_query_ratio-observer.update_query_ratio>=MIN_UQ_DIFF)&&(observer.update_query_ratio<=Update_Query_Threshold)){
                            if(tp[0].num_thread_update==big_u_update&&tp[1].num_thread_update==big_u_update){
                                mode = RESET_MODE;
                                r_query_num = big_q_query;
                                r_update_num = big_q_update;
                            }else if(tp[0].num_thread_update==big_u_update||tp[1].num_thread_update==big_u_update)
                            {
                                mode = EVALUATION_START_MODE;
                                last_eva_t = current_time;
                                start_evaluation();
                                cout<<"Last mode: NORMAL--> Start Evaluation!"<<endl<<endl;
                            }
                        }
                        if(mode == RESET_MODE)
                        {
                            restart_flag = 1;
                            update_query_time();
                            clear_query_time();
                            tp[0].restart_flag = 1;
                            cout << "Mode Last: NORMAL" << endl << "Reset Pool 0 To:" << r_query_num << "," << r_update_num
                                 << "!!" << endl << endl << endl;
                            task_reinit_withQU(0, r_query_num, r_update_num);
                            if (can_estimate)
                                gettimeofday(&end, NULL);
                            else {
                                estimate_mutex.lock();
                                gettimeofday(&end, NULL);
                                estimate_mutex.unlock();
                            }
                            current_time = (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec -
                                           global_start.tv_usec;
                            mode = EVALUATION_START_MODE;
                            last_eva_t = current_time;
                            start_evaluation();
//                            No_Record_Flag = 1;
                            Start_No_Record_Time = current_time;
                        }
                    }
                }
                else if(mode == EVALUATION_START_MODE){
                    if(current_time-last_eva_t>=EVA_TIME*MICROSEC_PER_SEC)
                    {
                        cout<<"End Evaluation!"<<endl;
                        mode = EVALUATION_END_MODE;
                        end_evaluation();
                        update_query_time();
                        clear_query_time();
                        int restart_pool = 0;
                        int go_to_flag = -1;
                        if(tp[0].eva_response_speed<tp[1].eva_response_speed)
                        {
                            restart_pool = 0;
                        }
                        else {
                            restart_pool = 1;
                        }
                        if(tp[1-restart_pool].num_thread_update==big_q_update) go_to_flag = QUERY_SET;
                        else go_to_flag = UPDATE_SET;
                        if(go_to_flag==QUERY_SET&&observer.update_query_ratio>=Update_Query_Threshold){
                            cout<<"change UQ-Threshold from:"<<Update_Query_Threshold<<" to "<<observer.update_query_ratio<<endl;
                            Update_Query_Threshold = observer.update_query_ratio;
                        }else if(go_to_flag==UPDATE_SET&&observer.update_query_ratio<=Update_Query_Threshold){
                            cout<<"change UQ-Threshold from:"<<Update_Query_Threshold<<" to "<<observer.update_query_ratio<<endl;
                            Update_Query_Threshold = observer.update_query_ratio;
                        }
                        tp[restart_pool].restart_flag = 1;
                        cout<<"restart pool"<<restart_pool<<"!!"<<endl<<endl<<endl;
                        task_reinit(restart_pool);
                        restart_flag = 1;
                        mode = NORMAL_MODE;
                        if (can_estimate)
                            gettimeofday(&end, NULL);
                        else {
                            estimate_mutex.lock();
                            gettimeofday(&end, NULL);
                            estimate_mutex.unlock();
                        }
                        current_time = (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec -
                                       global_start.tv_usec;
//                        No_Record_Flag = 1;
                        Start_No_Record_Time = current_time;
                    }
                }
                if(X_STAR_MODE)   compute_x_star();
                t_min+=1;
            }

//            if(No_Record_Flag &&(current_time-Start_No_Record_Time>=MAX_NORECORD_TIME*MICROSEC_PER_SEC))
//            {
//                No_Record_Flag = 0;
//            }



            if(i==0) offset_time = current_time-issue_time;
            if(restart_flag){
                offset_time = current_time-issue_time;
                restart_flag = 0;
            }
            if(!VERIFY) {

                do {
                    if (issue_time <= current_time-offset_time) {
                        break;
                    }
                    if(can_estimate)
                        gettimeofday(&end, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&end, NULL);
                        estimate_mutex.unlock();
                    }
                    current_time =
                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                    if (issue_time > current_time-offset_time) std::this_thread::sleep_for(std::chrono::microseconds(1));
                } while (true);
                issue_time = current_time;
            }

            if(observer.task_list.size()>=NUM_OBV_T*multiTestPara.init_objects)
            {
                observer.task_list.erase(observer.task_list.begin(),observer.task_list.begin()+1);
                observer.time_list.erase(observer.time_list.begin(),observer.time_list.begin()+1);
            }
            if(event.second==QUERY) {
                observer.task_list.push_back(QUERY);
            }
            else observer.task_list.push_back(1-QUERY);
            observer.time_list.push_back(current_time);

//            cout<<"step1: event-"<<event.second<<endl;
            // if insert
            if (event.second == INSERT) {

                int non_object_node = arrival_task_nodes[i];

//                if(need_opt){
//                    non_object_node%=1270000;
//                }
                for(int ti = 0;ti<2;ti++)
                {
//                    cout<<"id:"<<ti<<endl;
                    int start_q_id = tp[ti].global_start_q_id;
                    tp[ti].global_start_q_id=(tp[ti].global_start_q_id+1)%tp[ti].num_thread_update;
                    tp[ti].rand_idx_update=(tp[ti].rand_idx_update+1)%rand_length+rand_length;

                    for(int z = 0; z<tp[ti].num_threads_query;z++) {


                        pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                        pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);


                        // assign insert tasks
                        int min_task_size = INT_MAX;
                        int pool_index=z * tp[ti].num_thread_update + start_q_id;
                        tp[ti]._pool[pool_index]->add_task(task);
                        tp[ti].total_object_map[z][non_object_node] = pool_index;
                    }
                    tp[ti].current_object_node.push_back(non_object_node);
                }
                if(VERIFY){
//                    car_nodes[non_object_node]=1;
                    DijkstraKNNInsert(non_object_node, car_nodes);
                }
            }
            if (event.second == DELETE) {

                int object_node = arrival_task_nodes[i];
//                if(need_opt){
//                    object_node%=1270000;
//                }
                for(int ti=0;ti<2;ti++)
                {
                    for(int z = 0; z<tp[ti].num_threads_query;z++) {
                        int start_q_id = tp[ti].total_object_map[z][object_node];
                        for (int j = start_q_id; j < 1 + start_q_id; j++) {
                            int j_mod = j % tp[ti].num_thread_update;
                            pair<int, int> node_type_pair = std::make_pair(object_node, DELETE);
                            pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
                            tp[ti]._pool[z * tp[ti].num_thread_update + j_mod]->add_task(task);
                        }
                    }
                    int key_obj = -1;
                    for(int obj_i = 0;obj_i<tp[ti].current_object_node.size();obj_i++){
                        if(tp[ti].current_object_node[obj_i]==object_node){
                            key_obj = obj_i;
                            break;
                        }
                    }
                    if(key_obj!=-1){
                        int last_index = tp[ti].current_object_node.size() - 1;
                        int tmp = tp[ti].current_object_node[key_obj];
                        tp[ti].current_object_node[key_obj] = tp[ti].current_object_node[last_index];
                        tp[ti].current_object_node[last_index] = tmp;
                        tp[ti].current_object_node.pop_back();
                    }
                    else{
                        cout<<"not exist object node:"<<object_node<<endl<<endl<<endl;
                    }
                }
                if(VERIFY){
//                    car_nodes[object_node]=0;
                    DijkstraKNNDelete(object_node, car_nodes);
                }
            }
            if (event.second == QUERY) {
//                for(int qid =0;qid<2;qid++)
//                {
//                    tp[qid].ontask_num = 0;
//                    for(int z = 0;z < tp[qid].num_threads_query;z++)
//                    {
//                        for(int q_id = 0;q_id <tp[qid].num_thread_update;q_id++)
//                        {
//                            int pool_index=z * tp[qid].num_thread_update + q_id;
//                            int num_queries = tp[qid]._pool[pool_index]->get_num_queries_in_queue();
//                            int num_inserts = tp[qid]._pool[pool_index]->get_num_inserts_in_queue();
//                            int num_deletes = tp[qid]._pool[pool_index]->get_num_deletes_in_queue();
//                            tp[qid].ontask_num+= num_queries+num_inserts+num_deletes;
//                        }
//                    }
//                }
//                if(tp[0].ontask_num>tp[1].ontask_num)
//                {
//                    query_turn_flag = 1;
//                }else if(tp[0].ontask_num<tp[1].ontask_num) {
//                    query_turn_flag = 0;
//                }else{
//                    query_turn_flag = 1-query_turn_flag;
//                }
                query_turn_flag = 1-query_turn_flag;
//                int ti = query_turn_flag;
//                turn_num ++;
//                task_turn_mutex.lock();
//                if(turn_num%100==0)
//                {
//                    query_turn_flag = 1-query_turn_flag;
//                    cout<<"pool "<<query_turn_flag<<endl;
//                }
//                task_turn_mutex.unlock();
                int ti = query_turn_flag;
//                cout<<"pool "<<query_turn_flag<<endl;
//                int ti = 0;
                tp[ti].total_queries++;
                if(can_estimate) {
                    gettimeofday(&end, NULL);
                }
                else{
                    estimate_mutex.lock();
                    gettimeofday(&end, NULL);
                    estimate_mutex.unlock();
                }
                long current_time2 =
                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                tp[ti].total_offset+=current_time2-issue_time;

                int query_node=arrival_task_nodes[i];
//                if(need_opt){
//                    query_node%=1270000;
//                }

//                if(VERIFY){
//                    vector<KNode> result = DijkstraKNNQuery(k, query_node, mems.dist,
//                                                            mems.visited, mems.q, car_nodes);
//                    verify_results.push_back(result);
//
//                }
                // put to query tasks
                for (int j = 0; j < tp[ti].num_thread_update; j++) {
                    pair<int, int> node_type_pair = std::make_pair(query_node, QUERY);
                    pair<long, pair<int, int> > task = std::make_pair(current_time2, node_type_pair);
//                    cout<<"query added to "<<current_query_threads * num_threads_update + j<<endl;
                    if(multiTestPara.method_name.compare("dijk")!=0) {
                        int use_thread_id=-1;
                        long min_cost = INT_MAX;
                        for(int u = tp[ti].current_query_threads;u < tp[ti].current_query_threads+tp[ti].num_threads_query; u++) {
                            int u_mod = u % tp[ti].num_threads_query;
                            long est_cost = tp[ti]._pool[u_mod * tp[ti].num_thread_update + j]->get_est_cost();
                            if(min_cost>est_cost){
                                min_cost = est_cost;
                                use_thread_id = u_mod;
                            }
                        }
                        tp[ti]._pool[use_thread_id * tp[ti].num_thread_update + j]->add_task(task);
                    }
                    else
                        tp[ti]._pool[tp[ti].current_query_threads * tp[ti].num_thread_update + j]->add_task(task);
                }
                tp[ti].current_query_threads=(tp[ti].current_query_threads+1)%tp[ti].num_threads_query;

                if(can_estimate) {
                    gettimeofday(&end, NULL);
                }
                else{
                    estimate_mutex.lock();
                    gettimeofday(&end, NULL);
                    estimate_mutex.unlock();
                }
                long current_time3 =
                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
                if(observer.ts[ti].size()>EXP_SIZE)
                {
                    observer.ts[ti].erase(observer.ts[ti].begin(),observer.ts[ti].begin()+1);
                }
                observer.ts[ti].push_back(current_time3-issue_time);
//                cout<<"query assign cost: "<<clock()-start_1<<endl;
            }
//           if (i%1000 == 0)
//           {
//               cout<<"i:"<<i<<endl;
//               for(int z = 0;z < num_threads_query;z++)
//               {
//                   for(int q_id = 0;q_id <num_threads_update;q_id++)
//                   {
//                       int pool_index=z * num_threads_update + q_id;
//                       int num_queries = _pool[pool_index]->get_num_queries_in_queue();
//                       int num_inserts = _pool[pool_index]->get_num_inserts_in_queue();
//                       int num_deletes = _pool[pool_index]->get_num_deletes_in_queue();
//                       cout<<"query:"<<z<<" update:"<<q_id<<" queries:"<<num_queries<<" inserts:"<<num_inserts<<" deletes:"<<num_deletes<<endl;
//                   }
//               }
//           }

        }
    }
    void task_reinit_withQU(int ti,int query_num,int update_num)
    {
        int tj = 1-ti;
        wait_for_finish(ti);
        set_stop(ti);
        join(ti);
        update_query_time();
        observer.task_list.clear();
        observer.time_list.clear();
        observer.ts[0].clear();
        observer.ta[0].clear();
        observer.tu_ex[0].clear();
        observer.tq_ex[0].clear();
        observer.ts[1].clear();
        observer.ta[1].clear();
        observer.tu_ex[1].clear();
        observer.tq_ex[1].clear();
        observer.update_query_ratio_set.clear();
        long total_response_time = tp[0].response_time+tp[1].response_time;
        number_of_queries = tp[0].query_num+tp[1].query_num;
//        cout<<"toal response time 0 :"<<tp[0].response_time<<" number of queries 0: "<<tp[0].query_num<<endl;
//        cout<<"toal response time 1 :"<<tp[1].response_time<<" number of queries 1: "<<tp[1].query_num<<endl;
//        cout << "expected response time: " << total_response_time / float(number_of_queries) << " seconds" << endl;
//        cout << "total_response_time: " << total_response_time << endl;
//        cout << "number_of_queries: " << number_of_queries << endl;
        tp[0].response_time_first = tp[0].response_time;
        tp[0].query_num_first = tp[0].query_num;
        tp[1].response_time_first = tp[1].response_time;
        tp[1].query_num_first = tp[1].query_num;
        tp[ti].threadpool_id=ti;
        tp[ti].num_threads_query = query_num;
        tp[ti].num_thread_update = update_num;
        tp[ti].run_time++;
        tp[ti]._needjoin = 0;
        tp[ti].threshold_number = 0;

        tp[ti].total_object_map.clear();
        tp[ti].init_arrival_nodes.clear();
        tp[ti].init_list.clear();
        tp[ti].current_object_node.clear();
        tp[ti].init_arrival_nodes.assign(tp[tj].current_object_node.begin(),tp[tj].current_object_node.end());
        for(int init_i = 0;init_i< tp[ti].init_arrival_nodes.size();init_i++)
        {
            tp[ti].init_list.push_back(make_pair(0.0, INSERT));
        }


        tp[ti].current_frame = 0;
        tp[ti].begin_frame = 0;

        tp[ti].global_start_q_id = 0;
        tp[ti].restart_flag = 0;
        tp[ti].rand_idx_query = 0;
        tp[ti].rand_idx_update = 0;
        tp[ti].restart_flag = 0;
        tp[ti]._pool.clear();
//        cout<<"init step1"<<endl;

        globalThreadVar[ti] = new GlobalThreadVar*[tp[ti].num_threads_query];

        for(int id =0;id<2;id++)
        {
            cout<<"pool: "<<id<<"update: "<<tp[id].num_thread_update<<"query: "<<tp[id].num_threads_query<<endl;
        }
//        cout<<"init step2"<<endl;
        int k_star = compute_k_star(k, tp[ti].num_thread_update, alpha, fail_p);
        for(int j =0;j<tp[ti].num_threads_query;j++) {
            globalThreadVar[ti][j] = new GlobalThreadVar();
            globalThreadVar[ti][j]->ran_threshold.clear();
//            cout<<"query:"<<j<<" time:"<<globalThreadVar[j]->total_query_time<<" num:"<<globalThreadVar[j]->number_of_queries<<endl;
            for (int i = 0; i <= QUERY_ID_FOLD; i++) globalThreadVar[ti][j]->ran_threshold.push_back(INT_MAX);
        }
//        cout<<"init step3"<<endl;
        if(!multiTestPara.is_single_aggregate)
            tp[ti]._aggregate_thread = new RandomAggregateThread* [tp[ti].num_threads_query];
        if(multiTestPara.is_single_aggregate){
            tp[ti]._single_aggregate_thread=new RandomAggregateThread(ti,0, k, tp[ti].num_thread_update);
        }
        if(DISPLAY)cout << "k_star: " << k_star << endl;
        for(int j =0;j<tp[ti].num_threads_query;j++){

            if(!multiTestPara.is_single_aggregate)
                tp[ti]._aggregate_thread[j] = new RandomAggregateThread(ti,j, k, tp[ti].num_thread_update);
            for(int i=0;i<tp[ti].num_thread_update;i++) {
                if(!multiTestPara.is_single_aggregate) {
                    RandomThread *t = new RandomThread(ti,j, i, k_star,query_cost,insert_cost,delete_cost, tp[ti]._aggregate_thread[j]);
                    tp[ti]._pool.push_back(t);
                }
                else{
                    RandomThread *t = new RandomThread(ti,j, i, k_star,query_cost,insert_cost,delete_cost, tp[ti]._single_aggregate_thread);
                    tp[ti]._pool.push_back(t);
                }
            }
        }
        clear_query_time();
        task_init(ti);
    }
    void task_reinit(int ti)
    {
        int tj = 1-ti;
        wait_for_finish(ti);
        set_stop(ti);
        join(ti);
        update_query_time();
        observer.task_list.clear();
        observer.time_list.clear();
        observer.ts[0].clear();
        observer.ta[0].clear();
        observer.tu_ex[0].clear();
        observer.tq_ex[0].clear();
        observer.ts[1].clear();
        observer.ta[1].clear();
        observer.tu_ex[1].clear();
        observer.tq_ex[1].clear();
        long total_response_time = tp[0].response_time+tp[1].response_time;
        number_of_queries = tp[0].query_num+tp[1].query_num;
        cout<<"toal response time 0 :"<<tp[0].response_time<<" number of queries 0: "<<tp[0].query_num<<endl;
        cout<<"toal response time 1 :"<<tp[1].response_time<<" number of queries 1: "<<tp[1].query_num<<endl;
        cout << "expected response time: " << total_response_time / float(number_of_queries) << " seconds" << endl;
        cout << "total_response_time: " << total_response_time << endl;
        cout << "number_of_queries: " << number_of_queries << endl;
        tp[0].response_time_first = tp[0].response_time;
        tp[0].query_num_first = tp[0].query_num;
        tp[1].response_time_first = tp[1].response_time;
        tp[1].query_num_first = tp[1].query_num;
        tp[ti].threadpool_id=ti;
        tp[ti].num_threads_query = tp[tj].num_threads_query;
        tp[ti].num_thread_update = tp[tj].num_thread_update;
        tp[ti].run_time++;
        tp[ti]._needjoin = 0;
        tp[ti].threshold_number = 0;

        tp[ti].total_object_map.clear();
        tp[ti].init_arrival_nodes.clear();
        tp[ti].init_list.clear();
        tp[ti].current_object_node.clear();
        tp[ti].init_arrival_nodes.assign(tp[tj].current_object_node.begin(),tp[tj].current_object_node.end());
        for(int init_i = 0;init_i< tp[ti].init_arrival_nodes.size();init_i++)
        {
            tp[ti].init_list.push_back(make_pair(0.0, INSERT));
        }


        tp[ti].current_frame = 0;
        tp[ti].begin_frame = 0;

        tp[ti].global_start_q_id = 0;
        tp[ti].restart_flag = 0;
        tp[ti].rand_idx_query = 0;
        tp[ti].rand_idx_update = 0;
        tp[ti].restart_flag = 0;
        tp[ti]._pool.clear();
//        cout<<"init step1"<<endl;

        globalThreadVar[ti] = new GlobalThreadVar*[tp[ti].num_threads_query];

        for(int id =0;id<2;id++)
        {
            cout<<"pool: "<<id<<"update: "<<tp[id].num_thread_update<<"query: "<<tp[id].num_threads_query<<endl;
        }
//        cout<<"init step2"<<endl;
        int k_star = compute_k_star(k, tp[ti].num_thread_update, alpha, fail_p);
        for(int j =0;j<tp[ti].num_threads_query;j++) {
            globalThreadVar[ti][j] = new GlobalThreadVar();
            globalThreadVar[ti][j]->ran_threshold.clear();
//            cout<<"query:"<<j<<" time:"<<globalThreadVar[j]->total_query_time<<" num:"<<globalThreadVar[j]->number_of_queries<<endl;
            for (int i = 0; i <= QUERY_ID_FOLD; i++) globalThreadVar[ti][j]->ran_threshold.push_back(INT_MAX);
        }
//        cout<<"init step3"<<endl;
        if(!multiTestPara.is_single_aggregate)
            tp[ti]._aggregate_thread = new RandomAggregateThread* [tp[ti].num_threads_query];
        if(multiTestPara.is_single_aggregate){
            tp[ti]._single_aggregate_thread=new RandomAggregateThread(ti,0, k, tp[ti].num_thread_update);
        }
        if(DISPLAY)cout << "k_star: " << k_star << endl;
        for(int j =0;j<tp[ti].num_threads_query;j++){

            if(!multiTestPara.is_single_aggregate)
                tp[ti]._aggregate_thread[j] = new RandomAggregateThread(ti,j, k, tp[ti].num_thread_update);
            for(int i=0;i<tp[ti].num_thread_update;i++) {
                if(!multiTestPara.is_single_aggregate) {
                    RandomThread *t = new RandomThread(ti,j, i, k_star,query_cost,insert_cost,delete_cost, tp[ti]._aggregate_thread[j]);
                    tp[ti]._pool.push_back(t);
                }
                else{
                    RandomThread *t = new RandomThread(ti,j, i, k_star,query_cost,insert_cost,delete_cost, tp[ti]._single_aggregate_thread);
                    tp[ti]._pool.push_back(t);
                }
            }
        }
        clear_query_time();
        task_init(ti);
    }
    void wait_for_finish(int ti)
    {
        tp[ti].num_intask = 0;
        while(1) {
            for (int z = 0; z < tp[ti].num_threads_query; z++) {
                for (int q_id = 0; q_id < tp[ti].num_thread_update; q_id++) {
                    int pool_index = z * tp[ti].num_thread_update + q_id;
                    int num_queries = tp[ti]._pool[pool_index]->get_num_queries_in_queue();
                    int num_inserts = tp[ti]._pool[pool_index]->get_num_inserts_in_queue();
                    int num_deletes = tp[ti]._pool[pool_index]->get_num_deletes_in_queue();
                    tp[ti].num_intask += num_deletes + num_inserts + num_queries;
//                    cout << "query:" << z << " update:" << q_id << " queries:" << num_queries << " inserts:"
//                         << num_inserts << " deletes:" << num_deletes << endl;
                }
            }
            if(tp[ti].num_intask == 0) break;
            else tp[ti].num_intask = 0;
        }
    }
    void set_stop(int ti){
        for (int i = 0; i < tp[ti]._pool.size(); i++) {
            tp[ti]._pool[i]->set_stop();
            tp[ti]._pool[i]->notify();

        }
        if(!multiTestPara.is_single_aggregate) {
            for (int j = 0; j < num_threads_query; j++) {
                tp[ti]._aggregate_thread[j]->set_stop();
                tp[ti]._aggregate_thread[j]->notify();

            }
        }
        if(multiTestPara.is_single_aggregate){
            tp[ti]._single_aggregate_thread->set_stop();
            tp[ti]._single_aggregate_thread->notify();
        }
    }
    void start_evaluation(){
        for(int id =0;id<2;id++) {
            tp[id].eva_response_time = 0;
            tp[id].eva_query_num = 0;
            for (int i = 0; i < tp[id].num_threads_query; i++) {
                tp[id].eva_response_time += globalThreadVar[id][i]->total_query_time;
                tp[id].eva_query_num += globalThreadVar[id][i]->number_of_queries;
            }
        }
    }
    void end_evaluation(){
        for(int id =0;id<2;id++) {
            double last_response_time = tp[id].eva_response_time;
            int last_query_num = tp[id].eva_query_num;
            tp[id].eva_response_time = 0;
            tp[id].eva_query_num = 0;
            for (int i = 0; i < tp[id].num_threads_query; i++) {
                tp[id].eva_response_time += globalThreadVar[id][i]->total_query_time;
                tp[id].eva_query_num += globalThreadVar[id][i]->number_of_queries;
            }
            tp[id].eva_response_time -= last_response_time;
            tp[id].eva_query_num -= last_query_num;
            if(tp[id].eva_response_time<=1e-7 || tp[id].eva_query_num<=1e-7)
                tp[id].eva_response_speed = 100;
            else
                tp[id].eva_response_speed = tp[id].eva_query_num/tp[id].eva_response_time;
            cout<<"pool "<<id<<"response_speed "<<tp[id].eva_response_speed<<"response time "<<tp[id].eva_response_time<<"query num "<<tp[id].eva_query_num<<"last response time"<<last_response_time<<endl;
        }
    }
    void update_query_time(){
        for(int id =0;id<2;id++) {
            for (int i = 0; i < tp[id].num_threads_query; i++) {
                tp[id].response_time += globalThreadVar[id][i]->total_query_time;
                tp[id].query_num += globalThreadVar[id][i]->number_of_queries;
            }
        }
    }
    void clear_query_time(){
        for(int id =0;id<2;id++) {
            for (int i = 0; i < tp[id].num_threads_query; i++) {
                globalThreadVar[id][i]->total_query_time=0;
                globalThreadVar[id][i]->number_of_queries=0;
            }
        }
    }
    void Generate_results(){
        long total_response_time = tp[0].response_time+tp[1].response_time;
        number_of_queries = tp[0].query_num+tp[1].query_num;
        cout<<"toal response time 0 :"<<tp[0].response_time<<" number of queries 0: "<<tp[0].query_num<<endl;
        cout<<"toal response time 1 :"<<tp[1].response_time<<" number of queries 1: "<<tp[1].query_num<<endl;
        cout << "expected response time: " << total_response_time / float(number_of_queries) << " seconds" << endl;
        cout << "total_response_time: " << total_response_time << endl;
        cout << "number_of_queries: " << number_of_queries << endl;
        cout << "expected_update_response_time: " << total_update_response_time / number_of_updates << endl;

        cout << "expected_update_process_time: " << total_update_process_time / number_of_updates << endl;
        std::ofstream outfile;


        outfile.open(input_parameters.output_data_dir + "stone9_24_outfile_auto_threadeach_"+to_string(num_threads_each) +"_"+ configstr+".txt", std::ios_base::app);
        outfile << endl
                <<" toal response time 0 first: "<<tp[0].response_time_first<<" number of queries 0: "<<tp[0].query_num_first
                <<" toal response time 1 first: "<<tp[1].response_time_first<<" number of queries 1: "<<tp[1].query_num_first
                <<" query response time first:"<<(tp[0].response_time_first+tp[1].response_time_first)/float(tp[0].query_num_first+tp[1].query_num_first)
                <<" toal response time 0 :"<<tp[0].response_time<<" number of queries 0: "<<tp[0].query_num
                <<" toal response time 1 :"<<tp[1].response_time<<" number of queries 1: "<<tp[1].query_num
//                << network_name<<" "
//                << "init: "<<multiTestPara.init_objects<<" "
//                << multiTestPara.method_name << " config simulation time: "
//                << multiTestPara.config_simulation_time << " test simulate time: "
//                << multiTestPara.test_simulation_time << " configure: "
//                << configurationId << " threshold: " << multiTestPara.is_thresholded << " fail_p: " << fail_p
//                << " "
                << "query "<<configstr
//                << " " << delete_rate << " " << multiTestPara.method_name << " singleAggregate: "
//                << multiTestPara.is_single_aggregate << " "
//                << multiTestPara.num_threads_update
//                << " " << multiTestPara.num_threads_query << " "
                << "query response time: " << total_response_time / float(number_of_queries) << " "
                << "query process time: " << total_query_process_time / number_of_query_processings << " "
                << "update response time: " << total_update_response_time / number_of_updates << " "
                << "update process time: " << total_update_process_time / number_of_updates
//                << " overload: " << overload_flag
                <<" query finish: "<<query_finish_rate
                <<" update finish: "<<1.0-update_finish_rate
                <<" schedule cost: "<<avg_offset
                << endl;
        outfile.close();
    }
    void run() {
        task_init(0);
        task_init(1);
        long offset_time;
        struct timeval global_start_2;
        if(can_estimate)
            gettimeofday(&global_start_2, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&global_start_2, NULL);
            estimate_mutex.unlock();
        }
        task_run();
        wait_for_finish(0);
        wait_for_finish(1);
        struct timeval end_2;
        if(can_estimate)
            gettimeofday(&end_2, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&end_2, NULL);
            estimate_mutex.unlock();
        }

        long duration =
                (end_2.tv_sec - global_start_2.tv_sec);
        if(DISPLAY){
            cout<<"duration: "<<duration<<" secs; fulllist size: "<<full_task_list.size()+tp[0].init_list.size()<<endl;
            avg_offset = tp[0].total_offset/(tp[0].total_queries+1);
            cout<<"Pool 0 avg offset: "<<avg_offset<<"queries: "<<tp[0].total_queries<<endl;
            avg_offset = tp[1].total_offset/(tp[1].total_queries+1);
            cout<<"Pool 1 avg offset: "<<avg_offset<<"queries: "<<tp[1].total_queries<<endl;
        }
        update_query_time();
        set_stop(0);
        set_stop(1);
        tp[0]._needjoin = 1;
        tp[1]._needjoin = 1;
        if(DISPLAY)cout << "all set stopped!" << endl;
        join_all();
        Generate_results();
    }
};
class RandomThreadPool_new {
private:
    int threadpool_id;
    int begin_node;
    int end_node;
    vector<RandomThread *> _pool;
    std::thread _main_thread;
    RandomAggregateThread **_aggregate_thread;

    RandomAggregateThread * _single_aggregate_thread;
    int num_threads_update;
    int _needjoin;
    double alpha;
    int test_n;
    double fail_p;
    int k;
    // stores the part where the object locates
    vector< vector<int> >total_object_map;
    vector<int> distribute_sizes;

    double query_rate;
    double insert_rate;
    double delete_rate;
    int num_threads_query;
    int simulation_time;
    long total_queries_plan;
    long total_updates_plan;
    long total_queries_finished;
    long total_updates_finished;
    vector<std::pair<double, int> > full_task_list;
    vector<std::pair<double, int> > init_list;
    vector<int> arrival_task_nodes;
    vector<int> init_arrival_nodes;
    int query_cost;
    int insert_cost;
    int delete_cost;
    int run_time; //how many times have run
    int threshold_number;
    vector<int> current_object_node;

    int i;//frame
    int begin_i;//begin_frame

    // Algorithm data structure
//    vector<int *> dijkstra_object_map_vec;
public:
    RandomThreadPool_new(int threadpool_id_val, int begin_node_val, int end_node_val, int num_threads_query_val, int num_threads_update_val, double alpha_val, int k_val, double fail_p_val,
                         int test_n_val, double query_rate_val, double insert_rate_val, double delete_rate_val, int simulation_time_val,int query_cost_val,int insert_cost_val,int delete_cost_val,
                         vector<std::pair<double, int> > full_task_list_val,vector<int> arrival_task_nodes_val,vector<std::pair<double, int> > init_list_val,vector<int> init_arrival_nodes_val,int run_time_val,
                         int threshold_number_val,int begin_i_val) {
        cout << "constructing RandomThreadPool..." << endl;
        if(overload_flag)
            overload_flag=0;
        i = begin_i_val;
        begin_i = begin_i_val;
        current_object_node.clear();
        threadpool_id = threadpool_id_val;
        begin_node = begin_node_val;
        end_node = end_node_val;
        simulation_time = simulation_time_val;
        num_threads_query=num_threads_query_val; //replication
        query_rate = query_rate_val;
        insert_rate = insert_rate_val;
        delete_rate = delete_rate_val;
        test_n = test_n_val;
        num_threads_update = num_threads_update_val; //partition
        if(DISPLAY)cout<<"num_threads_update:"<<num_threads_update<<" num_threads_query:"<<num_threads_query_val<<endl;
        alpha = alpha_val;
        fail_p = fail_p_val;
        k = k_val;
        total_queries_plan=0;
        total_updates_plan=0;
        total_queries_finished=0;
        total_updates_finished=0;
        full_task_list.assign(full_task_list_val.begin(),full_task_list_val.end());
        arrival_task_nodes.assign(arrival_task_nodes_val.begin(),arrival_task_nodes_val.end());
        init_list.assign(init_list_val.begin(),init_list_val.end());
        init_arrival_nodes.assign(init_arrival_nodes_val.begin(),init_arrival_nodes_val.end());
        query_cost = query_cost_val;
        insert_cost = insert_cost_val;
        delete_cost = delete_cost_val;
        _needjoin = 0;
        run_time = run_time_val;
        threshold_number = threshold_number_val;
        globalThreadVar[threadpool_id] = new GlobalThreadVar*[num_threads_query];
        int k_star = compute_k_star(k, num_threads_update, alpha, fail_p);
        for(int j =0;j<num_threads_query;j++) {
            globalThreadVar[threadpool_id][j] = new GlobalThreadVar();
            globalThreadVar[threadpool_id][j]->ran_threshold.clear();
//            cout<<"query:"<<j<<" time:"<<globalThreadVar[j]->total_query_time<<" num:"<<globalThreadVar[j]->number_of_queries<<endl;
            for (int i = 0; i <= QUERY_ID_FOLD; i++) globalThreadVar[threadpool_id][j]->ran_threshold.push_back(INT_MAX);
        }
        if(!multiTestPara.is_single_aggregate)
            _aggregate_thread = new RandomAggregateThread*[num_threads_query];
        if(multiTestPara.is_single_aggregate){
            _single_aggregate_thread=new RandomAggregateThread(threadpool_id,0, k, num_threads_update);
        }
        if(DISPLAY)cout << "k_star: " << k_star << endl;
        for(int j =0;j<num_threads_query;j++){

            if(!multiTestPara.is_single_aggregate)
                _aggregate_thread[j] = new RandomAggregateThread(threadpool_id,j, k, num_threads_update);
            for(int i=0;i<num_threads_update;i++) {
                if(!multiTestPara.is_single_aggregate) {
                    RandomThread *t = new RandomThread(threadpool_id,j, i, k_star,query_cost,insert_cost,delete_cost, _aggregate_thread[j]);
                    _pool.push_back(t);
                }
                else{
                    RandomThread *t = new RandomThread(threadpool_id,j, i, k_star,query_cost,insert_cost,delete_cost, _single_aggregate_thread);
                    _pool.push_back(t);
                }
            }
        }
    }

    //释放线程池
    ~RandomThreadPool_new() {
        if(DISPLAY)cout << "hi, RandomThreadPooladPool()" << endl;
//         for(int i = 0;i < _pool.size(); ++i){
//             delete _pool[i];
//         }

    }

    void join() {
        if(DISPLAY)cout << "start joining" << endl;

        if(!multiTestPara.is_single_aggregate) {
            for (int j = 0; j < num_threads_query; j++) {
                _aggregate_thread[j]->join();
            }
        }
        if(multiTestPara.is_single_aggregate){
            _single_aggregate_thread->join();
        }
        if(DISPLAY)cout << "finish joining aggregatethreads" << endl;



        int max_updates=0;
        for (int i = 0; i < _pool.size(); i++) {
            int remain_updates = _pool[i]->get_num_inserts_in_queue()+_pool[i]->get_num_deletes_in_queue();
            if(remain_updates > max_updates){
                max_updates =  remain_updates;
            }
            _pool[i]->join();
        }

        update_finish_rate = max_updates*1.0/total_updates_plan;
        for (int i = 0; i < num_threads_query; i++) {
            total_queries_finished += globalThreadVar[threadpool_id][i]->number_of_queries;
        }
        query_finish_rate = total_queries_finished * 1.0 / total_queries_plan;

        if(DISPLAY)cout << "finish joining threads" << endl;

        if(_main_thread.joinable())
            _main_thread.join();

        if(DISPLAY)cout << "finish joining main thread" << endl;

    }

    int isNeedJoin() {
        return _needjoin;
    }

    void start() {
        _main_thread = std::thread(&RandomThreadPool_new::run, this);
    }


    long get_max_remain_time(){
        if(isNeedJoin()) {
            long max_remain_time = 0;
            for (int i = 0; i < _pool.size(); i++) {
                long est_cost = _pool[i]->get_est_cost();
                if (est_cost > max_remain_time) {
                    max_remain_time = est_cost;
                }
            }
            return max_remain_time;
        } else
            return -1;

    }

    vector<int> get_current_object(){
            return current_object_node;
    };

    int get_current_tasknum(){
        return i;
    }
    void run() {
        mem_struct mems;
        int* car_nodes;
        if(VERIFY){
            allocate_mem(mems, test_n+1);
            car_nodes = new int[test_n+1];
            memset(car_nodes, 0, sizeof(int)*(test_n+1));
        }
        struct timeval start;
        struct timeval end;
        double alpha_floor = floor(alpha);
        double alpha_ceil = ceil(alpha);
        double r = rand() % RAND_MAX;
        double num_queues_selected = alpha_ceil;
        if (r < alpha - alpha_floor) {
            num_queues_selected = alpha_floor;
        }
        if(DISPLAY)cout << "num_queues_selected: " << num_queues_selected << endl;
        vector<int> object_list;
        vector<int> non_object_list;
        for(int j =0;j<40;j++) {
            vector<int> tmp;
            for (int i = 0; i <= test_n; i++) {
                tmp.push_back(0);
            }
            total_object_map.push_back(tmp);
        }
        for (int i = 0; i < test_n; i++) {
            non_object_list.push_back(i);

        }

        double last_time = 0.0;

        int current_query_threads=0;
        int rand_idx_query=0;
        int rand_idx_update=rand_length;
        int global_start_q_id=0;
        long total_offset=0;
        int total_queries=0;
        int need_opt=0;
        if(multiTestPara.method_name.compare("vtree")==0 && network_name.compare("BJ-old")==0){
            need_opt=1;
        }
//        int i;
//        cout<<"init:"<<init_objects<<endl;
        for(int init_i=0;init_i<init_list.size();init_i++){
//            cout<<init_i<<endl;
            if(can_estimate)
                gettimeofday(&end, NULL);
            else{
                estimate_mutex.lock();
                gettimeofday(&end, NULL);
                estimate_mutex.unlock();
            }
            long current_time =
                    (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

            pair<double, int> &event = init_list[init_i];
            int non_object_node = init_arrival_nodes[init_i];

            if(need_opt){
                non_object_node%=1270000;
            }

            int start_q_id = global_start_q_id;
            global_start_q_id=(global_start_q_id+1)%num_threads_update;
            rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;

            for(int z = 0; z<num_threads_query;z++) {


                pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                pair<long, pair<int, int> > task = std::make_pair(current_time, node_type_pair);


                // assign insert tasks
                int min_task_size = INT_MAX;

                int pool_index=z * num_threads_update + start_q_id;
                _pool[pool_index]->add_task(task);
                total_object_map[z][non_object_node] = pool_index;
            }
            current_object_node.push_back(non_object_node);
            if(VERIFY){
//                    car_nodes[non_object_node]=1;
                DijkstraKNNInsert(non_object_node, car_nodes);
            }

        }
//        if(run_time==1)
//        {
//            if(can_estimate)
//                gettimeofday(&global_start, NULL);
//            else{
//                estimate_mutex.lock();
//                gettimeofday(&global_start, NULL);
//                estimate_mutex.unlock();
//            }
//        }


        struct timeval global_start_2;
        if(can_estimate)
            gettimeofday(&global_start_2, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&global_start_2, NULL);
            estimate_mutex.unlock();
        }
        long offset_time;
        for (i=begin_i; i < full_task_list.size(); i++) {
            if(overload_flag) break;

            if(run_time ==0 && i>threshold_number)
            {
                _needjoin = 1;
                cout<<"Over threshold num!!!!!RESTART NOW!!!"<<endl<<endl<<endl<<endl<<endl;
                break;

            }
//            cout<<i<<endl;
            if(arrival_task_nodes[i]==-1) continue;
            pair<double, int> &event = full_task_list[i];
            long issue_time = floor(event.first * MICROSEC_PER_SEC);
//            if(event.second == QUERY){
//                if(threadpool_id!=task_turn_flag || i==current_task_num) continue;
//                else
//                {
//                    current_task_num = i;
//                    task_turn_flag = (task_turn_flag+1)%2;
//                }
//            }
//            cout<<"time: "<<event.first<<" sec"<<endl;
            if(can_estimate)
                gettimeofday(&end, NULL);
            else{
                estimate_mutex.lock();
                gettimeofday(&end, NULL);
                estimate_mutex.unlock();
            }
            long current_time=(end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
            if(i==begin_i) offset_time = current_time-issue_time;
            if(!VERIFY) {

                do {
                    if (issue_time <= current_time-offset_time) {
//                        if(event.second==QUERY && simulation_time==200) {
//                            if(issue_time + 50 < current_time) {
//                                cout << "current time > issue time + 50" << endl;
//                                cout << "current time: " << current_time << endl;
//                                cout << "issue time : " << issue_time << endl;
//                            }
//                        }
                        break;
                    }
                    if(can_estimate)
                        gettimeofday(&end, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&end, NULL);
                        estimate_mutex.unlock();
                    }
                    current_time =
                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                    if (issue_time > current_time) std::this_thread::sleep_for(std::chrono::microseconds(1));
                } while (true);
//                if(i>=begin_i){
                    issue_time = current_time;
//                }
//                if(i==begin_i&&run_time==0){
//                    if(can_estimate)
//                        gettimeofday(&global_start, NULL);
//                    else{
//                        estimate_mutex.lock();
//                        gettimeofday(&global_start, NULL);
//                        estimate_mutex.unlock();
//                    }
//                }
                // start from queue id $start_q_id, we list num_queues_selected consecutive queues to hold random updates
            }


            // if insert
            if (event.second == INSERT) {


                int non_object_node = arrival_task_nodes[i];

                if(need_opt){
                    non_object_node%=1270000;
                }



//                int start_q_id = global_random_numbers[rand_idx_update] % num_threads_update;
                int start_q_id = global_start_q_id;
                global_start_q_id=(global_start_q_id+1)%num_threads_update;
                rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;

                for(int z = 0; z<num_threads_query;z++) {


                    pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                    pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);


                    // assign insert tasks
                    int min_task_size = INT_MAX;
                    int pool_index=z * num_threads_update + start_q_id;
                    _pool[pool_index]->add_task(task);
                    total_object_map[z][non_object_node] = pool_index;
                }
                current_object_node.push_back(non_object_node);
                if(VERIFY){
//                    car_nodes[non_object_node]=1;
                    DijkstraKNNInsert(non_object_node, car_nodes);
                }
//                cout<<"insert assign cost : "<<clock()-start_1<<endl;
//                cout<<"end insert"<<endl;
            }
            if (event.second == DELETE) {

                int object_node = arrival_task_nodes[i];
                if(need_opt){
                    object_node%=1270000;
                }

                for(int z = 0; z<num_threads_query;z++) {
                    int start_q_id = total_object_map[z][object_node];
                    for (int j = start_q_id; j < 1 + start_q_id; j++) {
                        int j_mod = j % num_threads_update;
                        pair<int, int> node_type_pair = std::make_pair(object_node, DELETE);
                        pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
//                        cout<<"delete added to "<<z * num_threads_update + j_mod<<endl;
                        _pool[z * num_threads_update + j_mod]->add_task(task);
                    }
                }
                int key_obj = -1;
                for(int obj_i = 0;obj_i<current_object_node.size();obj_i++){
                    if(current_object_node[obj_i]==object_node){
                        key_obj = obj_i;
                        break;
                    }
                }
                if(key_obj!=-1){
                    int last_index = current_object_node.size() - 1;
                    int tmp = current_object_node[key_obj];
                    current_object_node[key_obj] = current_object_node[last_index];
                    current_object_node[last_index] = tmp;
                    current_object_node.pop_back();
                }
                else{
                    cout<<"not exist object node:"<<object_node<<endl<<endl<<endl;
                }
                if(VERIFY){
//                    car_nodes[object_node]=0;
                    DijkstraKNNDelete(object_node, car_nodes);
                }


            }
            if (event.second == QUERY) {

//                cout<<"query "<<endl;
                total_queries++;
                if(can_estimate) {
                    gettimeofday(&end, NULL);
                }
                else{
                    estimate_mutex.lock();
                    gettimeofday(&end, NULL);
                    estimate_mutex.unlock();
                }
                current_time =
                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                total_offset+=current_time-issue_time;

                int query_node=arrival_task_nodes[i];
                if(need_opt){
                    query_node%=1270000;
                }

                if(VERIFY){
                    vector<KNode> result = DijkstraKNNQuery(k, query_node, mems.dist,
                                                            mems.visited, mems.q, car_nodes);
                    verify_results.push_back(result);

                }
                // put to query tasks
                for (int j = 0; j < num_threads_update; j++) {
                    pair<int, int> node_type_pair = std::make_pair(query_node, QUERY);
                    pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
//                    cout<<"query added to "<<current_query_threads * num_threads_update + j<<endl;
                    if(multiTestPara.method_name.compare("dijk")!=0) {
                        int use_thread_id=-1;
                        long min_cost = INT_MAX;
                        for(int u = current_query_threads;u < current_query_threads+num_threads_query; u++) {
                            int u_mod = u % num_threads_query;
                            long est_cost = _pool[u_mod * num_threads_update + j]->get_est_cost();
                            if(min_cost>est_cost){
                                min_cost = est_cost;
                                use_thread_id = u_mod;
                            }

                        }

                        _pool[use_thread_id * num_threads_update + j]->add_task(task);
                    }
                    else
                        _pool[current_query_threads * num_threads_update + j]->add_task(task);

                }

                current_query_threads=(current_query_threads+1)%num_threads_query;
//                cout<<"query assign cost: "<<clock()-start_1<<endl;
            }
//           if (i%1000 == 0)
//           {
//               cout<<"i:"<<i<<endl;
//               for(int z = 0;z < num_threads_query;z++)
//               {
//                   for(int q_id = 0;q_id <num_threads_update;q_id++)
//                   {
//                       int pool_index=z * num_threads_update + q_id;
//                       int num_queries = _pool[pool_index]->get_num_queries_in_queue();
//                       int num_inserts = _pool[pool_index]->get_num_inserts_in_queue();
//                       int num_deletes = _pool[pool_index]->get_num_deletes_in_queue();
//                       cout<<"query:"<<z<<" update:"<<q_id<<" queries:"<<num_queries<<" inserts:"<<num_inserts<<" deletes:"<<num_deletes<<endl;
//                   }
//               }
//           }

        }
        int num_intask = 0;
        while(1) {
            for (int z = 0; z < num_threads_query; z++) {
                for (int q_id = 0; q_id < num_threads_update; q_id++) {
                    int pool_index = z * num_threads_update + q_id;
                    int num_queries = _pool[pool_index]->get_num_queries_in_queue();
                    int num_inserts = _pool[pool_index]->get_num_inserts_in_queue();
                    int num_deletes = _pool[pool_index]->get_num_deletes_in_queue();
                    num_intask += num_deletes + num_inserts + num_queries;
//                    cout << "query:" << z << " update:" << q_id << " queries:" << num_queries << " inserts:"
//                         << num_inserts << " deletes:" << num_deletes << endl;
                }
            }
            if(num_intask == 0) break;
            else num_intask = 0;
        }
        while(globalThreadVar[threadpool_id][0]->number_of_queries<2){
            std::this_thread::sleep_for(std::chrono::microseconds(1));
        }
        struct timeval end_2;
        if(can_estimate)
            gettimeofday(&end_2, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&end_2, NULL);
            estimate_mutex.unlock();
        }


        long duration =
                (end_2.tv_sec - global_start_2.tv_sec);
        if(DISPLAY)cout<<"duration: "<<duration<<" secs; fulllist size: "<<full_task_list.size()+init_list.size()<<endl;
        avg_offset = total_offset/total_queries;
        if(DISPLAY)cout<<"avg offset: "<<total_offset/total_queries<<endl;
        if(VERIFY){
            delete[] car_nodes;
            delete_mems(mems);
        }
//        for(int j = 0; j < num_threads_query;j++) {
        for (int i = 0; i < _pool.size(); i++) {
            _pool[i]->set_stop();
            _pool[i]->notify();

        }
        if(!multiTestPara.is_single_aggregate) {
            for (int j = 0; j < num_threads_query; j++) {
                _aggregate_thread[j]->set_stop();
                _aggregate_thread[j]->notify();

            }
        }
        if(multiTestPara.is_single_aggregate){
            _single_aggregate_thread->set_stop();
            _single_aggregate_thread->notify();
        }
        _needjoin = 1;
        if(DISPLAY)cout << "all set stopped!" << endl;
    }

};

class RandomThreadPool_Control{
private:
    int threadpool_id;
    int begin_node;
    int end_node;
    int num_threads;
    int num_threads_query[2];
    int num_threads_update[2];
    int _needjoin;
    double alpha;
    int test_n;
    double fail_p;
    int k;
    double query_rate;
    double insert_rate;
    double delete_rate;
    int simulation_time;
    long total_queries_plan;
    long total_updates_plan;
    long total_queries_finished;
    long total_updates_finished;
    int query_cost;
    int insert_cost;
    int delete_cost;
    struct timeval start, end;
    double total_response_time;
    RandomThreadPool_new * tp[2];
    int configurationId;
    int init_objects;
    vector<std::pair<double, int> > full_list;
    vector<int> arrival_nodes;
    vector<std::pair<double, int> > full_task_list;
    vector<int> arrival_task_nodes;
    int run_time[2]; //how many times we run x
    int threshold_number; //bigger than the threshold then we may restart
    int response_time[2]={0};
    int query_num[2] = {0};

public:
    RandomThreadPool_Control(int threadpool_id_val, int begin_node_val, int end_node_val, int num_threads_val, double alpha_val, int k_val, double fail_p_val,
                     int test_n_val, double query_rate_val, double insert_rate_val, double delete_rate_val, int simulation_time_val,int query_cost_val,int insert_cost_val,int delete_cost_val,int ConfigID_val){

        begin_node = begin_node_val;
        end_node = end_node_val;
        num_threads = num_threads_val;
        alpha = alpha_val;
        test_n = test_n_val;
        fail_p = fail_p_val;
        k = k_val;
        query_rate = query_rate_val;
        insert_rate = insert_rate_val;
        delete_rate = delete_rate_val;
        simulation_time = simulation_time_val;
        query_cost = query_cost_val;
        delete_cost = delete_cost_val;
        insert_cost = insert_cost_val;
        total_response_time = 0.0;
        configurationId = ConfigID_val;
        threshold_number = 400000;
}
    ~RandomThreadPool_Control(){

    };
    vector<int> generate_arrival_nodes(vector<std::pair<double, int> >& full_list, int begin, int end){

        int rand_idx_update=0;
        int rand_idx_query=0;
        vector<int> object_list;
        vector<int> non_object_list;
        vector<int> arrival_node_list;
        for (int i = begin; i < end; i++) {
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
                int index = rand() % non_object_list.size();
                rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;
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
                int index = rand() % object_list.size();
                rand_idx_update=(rand_idx_update+1)%rand_length;
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
                    rand_idx_query=(rand_idx_query+1)%rand_length;
                    query_node = query_node % (test_n-1)+1;

                } while (vStart[query_node] == 0);
                arrival_node_list.push_back(query_node);
            }

        }
        return arrival_node_list;
    }
    void Generate_results(){
//        if (!multiTestPara.is_single_aggregate) {
//            for (int i = 0; i < multiTestPara.num_threads_query; i++) {
//                total_response_time += globalThreadVar[i]->total_query_time;
//                number_of_queries += globalThreadVar[i]->number_of_queries;
//            }
//        } else {
//            total_response_time += globalThreadVar[0]->total_query_time;
//            number_of_queries += globalThreadVar[0]->number_of_queries;
//        }
        total_response_time = response_time[0]+response_time[1];
        number_of_queries = query_num[0]+query_num[1];
        cout<<"toal response time 0 :"<<response_time[0]<<" number of queries 0: "<<query_num[0]<<endl;
        cout<<"toal response time 1 :"<<response_time[1]<<" number of queries 1: "<<query_num[1]<<endl;
        cout << "expected response time: " << total_response_time / number_of_queries << " seconds" << endl;
        cout << "total_response_time: " << total_response_time << endl;
        cout << "number_of_queries: " << number_of_queries << endl;
        cout << "expected_update_response_time: " << total_update_response_time / number_of_updates << endl;

        cout << "expected_update_process_time: " << total_update_process_time / number_of_updates << endl;
        std::ofstream outfile;

        outfile.open(input_parameters.output_data_dir + "stone_outfile" + (multiTestPara.suffix), std::ios_base::app);
        outfile << endl
                << network_name<<" "
                << "init: "<<multiTestPara.init_objects<<" "
                << multiTestPara.method_name << " config simulation time: "
                << multiTestPara.config_simulation_time  << " configure: "
                << configurationId << " threshold: " << multiTestPara.is_thresholded << " fail_p: " << fail_p
                << " "
                << query_rate << " " << insert_rate
                << " " << delete_rate << " " << multiTestPara.method_name << " singleAggregate: "
                << multiTestPara.is_single_aggregate << " "
                << multiTestPara.num_threads_update
                << " " << multiTestPara.num_threads_query << " "
                << "query response time: " << total_response_time / number_of_queries << " "
                << "query process time: " << total_query_process_time / number_of_query_processings << " "
                << "update response time: " << total_update_response_time / number_of_updates << " "
                << "update process time: " << total_update_process_time / number_of_updates << " "
                << end.tv_sec - start.tv_sec
                << " overload: " << overload_flag
                <<" query finish: "<<query_finish_rate
                <<" update finish: "<<1.0-update_finish_rate
                <<" schedule cost: "<<avg_offset
                << endl;
        outfile.close();
    }

    void init(){
        full_task_list.clear();
        full_list.clear();
        arrival_nodes.clear();
        arrival_task_nodes.clear();
        init_objects = multiTestPara.init_objects;
        for (int i = 0; i < init_objects; i++) {
            full_list.push_back(make_pair(0.0, INSERT));
        }

        std::ifstream queryfile;
        queryfile.open(input_parameters.input_data_dir + "query_" + std::to_string(query_rate) + "_" +
                       std::to_string(insert_rate) + "_" +
                       std::to_string(delete_rate) + "_" + std::to_string(simulation_time) + ".txt",
                       std::ios_base::in);
        double f1;
        int f2;
        if (!queryfile.is_open()) {
            cout << "can't load queryfile!" << endl;
            vector<std::pair<double, int> > append_list = make_online_query_update_list(query_rate, insert_rate,
                                                                                        delete_rate,
                                                                                        simulation_time);
            std::ofstream queryfile_w;
            queryfile_w.open(input_parameters.input_data_dir + "query_" + std::to_string(query_rate) + "_" +
                             std::to_string(insert_rate) + "_" +
                             std::to_string(delete_rate) + "_" + std::to_string(simulation_time) + ".txt",
                             std::ios_base::out);
            for (pair<double, int> &item : append_list) {
                full_list.push_back(item);
                queryfile_w << item.first << " " << item.second << endl;
            }
            queryfile_w.close();
            cout << "write to queryfile!" << endl;
        } else {
            while (!queryfile.eof()) {
                queryfile >> f1 >> f2;
//                cout<<f1<<" "<<f2<<endl;
                full_list.push_back(make_pair(f1, f2));
            }
            queryfile.close();
            cout << "read from queryfile!" << endl;
        }

        arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);

        full_task_list.assign(full_list.begin()+init_objects,full_list.end());
        arrival_task_nodes.assign(arrival_nodes.begin()+init_objects,arrival_nodes.end());
        vector<std::pair<double, int> > init_list;
        init_list.assign(full_list.begin(),full_list.begin()+init_objects);
        vector<int> init_arrival_node;
        init_arrival_node.assign(arrival_nodes.begin(),arrival_nodes.begin()+init_objects);
        num_threads_update[0] = 1;
        num_threads_query[0] = num_threads - 2;
        int num_q = int(sqrt(num_threads-2));
        int num_p = int((num_threads-2)/num_q);
        num_threads_query[1] = max(num_q,num_p);
        num_threads_update[1] = num_p+num_q-num_threads_query[1];
        tp[0] = new RandomThreadPool_new(0, 0, end_node, num_threads_query[0],num_threads_update[0], alpha, k, fail_p,test_n,
                                        query_rate, insert_rate,delete_rate, simulation_time,
                                        query_cost,insert_cost,delete_cost,full_task_list,arrival_task_nodes,init_list,init_arrival_node,
                                        run_time[0],threshold_number,0);
        tp[1] = new RandomThreadPool_new(1, 0, end_node, num_threads_query[1],num_threads_update[1], alpha, k, fail_p,test_n,
                                         query_rate, insert_rate,delete_rate, simulation_time,
                                         query_cost,insert_cost,delete_cost,full_task_list,arrival_task_nodes,init_list,init_arrival_node,
                                         run_time[1],threshold_number,0);
        cout << "full_list made..." << endl;
        cout << "full list size: " << full_list.size() << endl;
    }
    void re_init(int id){
        run_time[id]++;
        int begin_frame = tp[id]->get_current_tasknum();
        vector<int> init_object_nodes = tp[id]->get_current_object();
        vector<std::pair<double, int> > init_list;
        for(int init_i = 0;init_i< init_object_nodes.size();init_i++)
        {
            init_list.push_back(make_pair(0.0, INSERT));
        }
//        delete tp_x;
        int num_q = int(sqrt(num_threads-2));
        int num_p = int((num_threads-2)/num_q);
        num_threads_query[id] = max(num_q,num_p);
        num_threads_update[id] = num_p+num_q-num_threads_query[id];
//        num_threads_query_x = num_threads-2;
//        num_threads_update_x = 1;
        tp[id] = new RandomThreadPool_new(id, 0, end_node, num_threads_query[id],num_threads_update[id], alpha, k, fail_p,test_n,
                                        query_rate, insert_rate,delete_rate, simulation_time,
                                        query_cost,insert_cost,delete_cost,full_task_list,arrival_task_nodes,init_list,init_object_nodes,
                                        run_time[id],threshold_number,begin_frame);
    }
    void update_query_time(){
        for(int id =0;id<2;id++) {
            for (int i = 0; i < num_threads_query[id]; i++) {
                response_time[id] += globalThreadVar[id][i]->total_query_time;
                query_num[id] += globalThreadVar[id][i]->number_of_queries;
            }
        }
    }
    void run()
    {
        init();
//        tp_x = new RandomThreadPool_new(0, 0, end_node, num_threads_query,num_threads_update, alpha, k, fail_p,test_n,
//                                  query_rate, insert_rate,delete_rate, simulation_time,
//                                  query_cost,insert_cost,delete_cost,full_list,arrival_nodes,init_objects,
//                                  x_time,threshold_number);
        tp[0]->start();
//        tp[1]->start();
        if(can_estimate)
            gettimeofday(&global_start, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&global_start, NULL);
            estimate_mutex.unlock();
        }
        while (true) {
            std::this_thread::sleep_for(std::chrono::microseconds(1));
            if (tp[0]->isNeedJoin()) {
                tp[0]->join();
                break;
            }
//            if (tp[1]->isNeedJoin()) {
//                tp[1]->join();
//            }
//            if(tp[0]->isNeedJoin()&&tp[1]->isNeedJoin()) break;
        }
        update_query_time();

        //

        re_init(0);
        tp[0]->start();
//        if(can_estimate)
//            gettimeofday(&global_start, NULL);
//        else{
//            estimate_mutex.lock();
//            gettimeofday(&global_start, NULL);
//            estimate_mutex.unlock();
//        }
        while (true) {
            std::this_thread::sleep_for(std::chrono::microseconds(1));
            if (tp[0]->isNeedJoin()) {
                if(can_estimate)
                    gettimeofday(&start, NULL);
                else{
                    estimate_mutex.lock();
                    gettimeofday(&start, NULL);
                    estimate_mutex.unlock();
                }
                tp[0]->join();
                break;
            }
        }
        update_query_time();
        delete tp[0];


        //
        if(can_estimate)
            gettimeofday(&end, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&end, NULL);
            estimate_mutex.unlock();
        }
        cout << end.tv_sec << " " << global_start.tv_sec << endl;
        cout << "finish in : " << end.tv_sec - global_start.tv_sec << " secs" << endl;
        Generate_results();
    }
};


class RandomThreadPool {
private:
    int threadpool_id;
    int begin_node;
    int end_node;
    vector<RandomThread *> _pool;
    // std::recursive_mutex _locker;
    std::thread _main_thread;
    RandomAggregateThread **_aggregate_thread;

    RandomAggregateThread * _single_aggregate_thread;
    int num_threads_update;
    int _needjoin;
    double alpha;
    int test_n;
    double fail_p;
    int k;
    // stores the part where the object locates
    vector< vector<int> >total_object_map;
    vector<int> distribute_sizes;

    double query_rate;
    double insert_rate;
    double delete_rate;
    int num_threads_query;
    int simulation_time;
    long total_queries_plan;
    long total_updates_plan;
    long total_queries_finished;
    long total_updates_finished;
    string configstr;

    // Algorithm data structure
//    vector<int *> dijkstra_object_map_vec;
public:
    RandomThreadPool(int threadpool_id_val, int begin_node_val, int end_node_val, int num_threads_val, double alpha_val, int k_val, double fail_p_val,
                     int test_n_val, string configstr_val,int query_cost,int insert_cost,int delete_cost, int configurationId_val) {
        cout << "constructing RandomThreadPool..." << endl;
        if(overload_flag)
            overload_flag=0;
        threadpool_id = threadpool_id_val;
        configstr = configstr_val;
        begin_node = begin_node_val;
        end_node = end_node_val;

        test_n = test_n_val;
        num_threads_update = multiTestPara.update_thread; //partition
        num_threads_query=floor((num_threads_val-2)/num_threads_update); //replication
        alpha = alpha_val;
        fail_p = fail_p_val;
        k = k_val;
        total_queries_plan=0;
        total_updates_plan=0;
        total_queries_finished=0;
        total_updates_finished=0;
        init();
        cout<<"num_threads_query: "<<num_threads_query<<"num_update_query: "<<num_threads_update<<endl;
        globalThreadVar[threadpool_id] = new GlobalThreadVar*[num_threads_query];
        int k_star = compute_k_star(k, num_threads_update, alpha, fail_p);
        for(int j =0;j<num_threads_query;j++) {
            globalThreadVar[threadpool_id][j] = new GlobalThreadVar();
            globalThreadVar[threadpool_id][j]->ran_threshold.clear();
            for (int i = 0; i <= QUERY_ID_FOLD; i++) globalThreadVar[threadpool_id][j]->ran_threshold.push_back(INT_MAX);
        }
        if(!multiTestPara.is_single_aggregate)
            _aggregate_thread = new RandomAggregateThread*[num_threads_query];
        if(multiTestPara.is_single_aggregate){
            _single_aggregate_thread=new RandomAggregateThread(threadpool_id,0, k, num_threads_update);
        }
        cout << "k_star: " << k_star << endl;
        for(int j =0;j<num_threads_query;j++){

            if(!multiTestPara.is_single_aggregate)
                _aggregate_thread[j] = new RandomAggregateThread(threadpool_id,j, k, num_threads_update);
            for(int i=0;i<num_threads_update;i++) {
                if(!multiTestPara.is_single_aggregate) {
                    RandomThread *t = new RandomThread(threadpool_id,j, i, k_star,query_cost,insert_cost,delete_cost, _aggregate_thread[j]);
                    _pool.push_back(t);
                }
                else{
                    RandomThread *t = new RandomThread(threadpool_id,j, i, k_star,query_cost,insert_cost,delete_cost, _single_aggregate_thread);
                    _pool.push_back(t);
                }
            }
        }

        _needjoin = 0;
    }

    //释放线程池
    ~RandomThreadPool()  {
        if(DISPLAY)cout << "hi, RandomThreadPooladPool()" << endl;
        // releasePool();
        // for(int i = 0;i < _pool.size(); ++i){
        //     delete _pool[i];
        // }

    }

    void init() {
        if(DISPLAY) cout << "start init() function ..." << endl;

    }

    void join() {
        cout << "start joining" << endl;

        if(!multiTestPara.is_single_aggregate) {
            for (int j = 0; j < num_threads_query; j++) {
                _aggregate_thread[j]->join();
            }
        }
        cout << "start joining aggregatethreads" << endl;
        if(multiTestPara.is_single_aggregate){
            _single_aggregate_thread->join();
        }
        cout << "finish joining aggregatethreads" << endl;



        int max_updates=0;
        for (int i = 0; i < _pool.size(); i++) {
            int remain_updates = _pool[i]->get_num_inserts_in_queue()+_pool[i]->get_num_deletes_in_queue();
            if(remain_updates > max_updates){
                max_updates =  remain_updates;
            }
            _pool[i]->join();
        }

        update_finish_rate = max_updates*1.0/total_updates_plan;
        for (int i = 0; i < num_threads_query; i++) {
            total_queries_finished += globalThreadVar[threadpool_id][i]->number_of_queries;
        }
        query_finish_rate = total_queries_finished * 1.0 / total_queries_plan;

        cout << "finish joining threads" << endl;

//        if(_main_thread.joinable())
//            _main_thread.join();
//
//        cout << "finish joining main thread" << endl;

    }

    int isNeedJoin() {
        return _needjoin;
    }
    vector<int> generate_arrival_nodes(vector<std::pair<double, int> >& full_list, int begin, int end){
//        mem_struct mems;
//        int *dijkstra_object_map;
//        allocate_mem(mems, test_n + test_n + 1);
//        dijkstra_object_map = new int[test_n + 1];
//        memset(dijkstra_object_map, 0, sizeof(int) * (test_n + 1));
//        vector<KNode>* hier_local_knn_arr_single;
//        if(multiTestPara.method_name.compare("toain")==0){
//            hier_local_knn_arr_single= new vector<KNode>[test_n+1];
//        }
//        int query_id = 0;
//
        int rand_idx_update=0;
        int rand_idx_query=0;
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
                rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;
                int non_object_node = non_object_list[index];
                int last_index = non_object_list.size() - 1;

                int tmp = non_object_list[index];
                non_object_list[index] = non_object_list[last_index];
                non_object_list[last_index] = tmp;
                non_object_list.pop_back();

                object_list.push_back(non_object_node);
                arrival_node_list.push_back(non_object_node);
//                HandleInsert(0, 0, multiTestPara.method_name, k, non_object_node, dijkstra_object_map,
//                             mems,
//                             hier_local_knn_arr_single);

            }
            if (event.second == DELETE && object_list.size() > 0) {
                total_updates_plan++;
//                int index = global_random_numbers[rand_idx_update] % object_list.size();
                int index = rand() % object_list.size();
                rand_idx_update=(rand_idx_update+1)%rand_length;
                int object_node = object_list[index];
                int last_index = object_list.size() - 1;
                int tmp = object_list[index];
                object_list[index] = object_list[last_index];
                object_list[last_index] = tmp;
                object_list.pop_back();
                non_object_list.push_back(object_node);
                arrival_node_list.push_back(object_node);
//                HandleDelete(0, 0, multiTestPara.method_name, k, object_node, dijkstra_object_map,
//                             mems, hier_local_knn_arr_single);

            }
            if (event.second == QUERY) {
                total_queries_plan++;
                int query_node;
                do {
//                    query_node = global_random_numbers[rand_idx_query];
                    query_node = rand();
                    rand_idx_query=(rand_idx_query+1)%rand_length;
                    query_node = query_node % (test_n-1)+1;

                } while (vStart[query_node] == 0);
                arrival_node_list.push_back(query_node);
//                HandleQuery(0, 0, multiTestPara.method_name, k, query_node, mems.dist,
//                            mems.visited, mems.q,
//                            dijkstra_object_map, globalThreadVar[0]->ran_threshold, query_id,
//                            hier_local_knn_arr_single);
//                query_id = (query_id+1)%QUERY_ID_FOLD;

            }

        }
//        stop_here();
        return arrival_node_list;

    }
    void releasePool() {
//        for(int j=0;j<num_threads_query;j++) {
        for (int i = 0; i < _pool.size(); ++i) {
            delete _pool[i];
//            }
        }
    }

//    void start() {
//
//        _main_thread = std::thread(&RandomThreadPool::run, this);
//    }


    long get_max_remain_time(){
        if(isNeedJoin()) {
            long max_remain_time = 0;
            for (int i = 0; i < _pool.size(); i++) {
                long est_cost = _pool[i]->get_est_cost();
                if (est_cost > max_remain_time) {
                    max_remain_time = est_cost;
                }
            }
            return max_remain_time;
        } else
            return -1;

    }

    void Generate_results(){
        double total_response_time=0;
        double number_of_queries=0;
        cout<<"ok"<<endl;
        for (int i = 0; i < num_threads_query; i++) {
            cout<<i<<endl;
            total_response_time += globalThreadVar[threadpool_id][i]->total_query_time;
            number_of_queries += globalThreadVar[threadpool_id][i]->number_of_queries;
        }
        cout << "expected response time: " << total_response_time / float(number_of_queries) << " seconds" << endl;
        cout << "total_response_time: " << total_response_time << endl;
        cout << "number_of_queries: " << number_of_queries << endl;
        cout << "expected_update_response_time: " << total_update_response_time / number_of_updates << endl;

        cout << "expected_update_process_time: " << total_update_process_time / number_of_updates << endl;
        std::ofstream outfile;


        outfile.open(input_parameters.output_data_dir + "stone9_24_outfile_auto_query_"+to_string(num_threads_query)+"_update_"+to_string(num_threads_update) +"_"+ configstr+".txt", std::ios_base::app);
        outfile << endl
                << " update thread "<<num_threads_update
                << " query thread "<<num_threads_query
                //                << network_name<<" "
                //                << "init: "<<multiTestPara.init_objects<<" "
                //                << multiTestPara.method_name << " config simulation time: "
                //                << multiTestPara.config_simulation_time << " test simulate time: "
                //                << multiTestPara.test_simulation_time << " configure: "
                //                << configurationId << " threshold: " << multiTestPara.is_thresholded << " fail_p: " << fail_p
                //                << " "
                << "query "<<configstr
                //                << " " << delete_rate << " " << multiTestPara.method_name << " singleAggregate: "
                //                << multiTestPara.is_single_aggregate << " "
                //                << multiTestPara.num_threads_update
                //                << " " << multiTestPara.num_threads_query << " "
                << "query response time: " << total_response_time / float(number_of_queries) << " "
                << "query process time: " << total_query_process_time / number_of_query_processings << " "
                << "update response time: " << total_update_response_time / number_of_updates << " "
                << "update process time: " << total_update_process_time / number_of_updates
                //                << " overload: " << overload_flag
                <<" query finish: "<<query_finish_rate
                <<" update finish: "<<1.0-update_finish_rate
                <<" schedule cost: "<<avg_offset
                << endl;
        outfile.close();
        _needjoin = 2;
    }
    void wait_for_finish(void)
    {
        int num_intask = 0;
        while(1) {
            for (int z = 0; z < num_threads_query; z++) {
                for (int q_id = 0; q_id < num_threads_update; q_id++) {
                    int pool_index = z * num_threads_update + q_id;
                    int num_queries = _pool[pool_index]->get_num_queries_in_queue();
                    int num_inserts = _pool[pool_index]->get_num_inserts_in_queue();
                    int num_deletes = _pool[pool_index]->get_num_deletes_in_queue();
                    num_intask += num_deletes + num_inserts + num_queries;
//                    cout << "query:" << z << " update:" << q_id << " queries:" << num_queries << " inserts:"
//                         << num_inserts << " deletes:" << num_deletes << endl;
                }
            }
            if(num_intask == 0) break;
            else num_intask = 0;
        }
    }
    void run() {
        mem_struct mems;
        int* car_nodes;
        if(VERIFY){
            allocate_mem(mems, test_n+1);
            car_nodes = new int[test_n+1];
            memset(car_nodes, 0, sizeof(int)*(test_n+1));
        }
        struct timeval start;
        struct timeval end;
        double alpha_floor = floor(alpha);
        double alpha_ceil = ceil(alpha);
        double r = rand() % RAND_MAX;
        double num_queues_selected = alpha_ceil;
        if (r < alpha - alpha_floor) {
            num_queues_selected = alpha_floor;
        }
        cout << "num_queues_selected: " << num_queues_selected << endl;
        vector<int> object_list;
        vector<int> non_object_list;
        for(int j =0;j<40;j++) {
            vector<int> tmp;
            for (int i = 0; i <= test_n; i++) {
                tmp.push_back(0);
//                distribute_sizes.push_back(0);
            }
            total_object_map.push_back(tmp);
        }
        for (int i = 0; i < test_n; i++) {
            non_object_list.push_back(i);

        }
        int init_objects = multiTestPara.init_objects;
        vector<std::pair<double, int> > full_list;
        for (int i = 0; i < init_objects; i++) {
            full_list.push_back(make_pair(0.0, INSERT));
        }

        string myversion = "v1_"+std::to_string(multiTestPara.init_objects)+"_";
        std::ifstream queryfile;

        string queryfile_name = input_parameters.input_data_dir + "query_"+myversion + configstr+".txt";
        cout<<"query file name:"<<queryfile_name<<endl;
        queryfile.open(queryfile_name,std::ios_base::in);
        double f1;
        int f2;
        if (!queryfile.is_open()) {
            cout << "can't load queryfile!" << endl;
            vector<std::pair<double, int> > append_list = make_online_query_update_list_str(configstr,multiTestPara.init_objects,multiTestPara.layer);
            cout<<"Generated!"<<endl;
            std::ofstream queryfile_w;
            queryfile_w.open(queryfile_name,std::ios_base::out);
            for (pair<double, int> &item : append_list) {
                full_list.push_back(item);
                queryfile_w << item.first << " " << item.second << endl;
            }
            queryfile_w.close();
            cout << "write to queryfile!" << endl;
        } else {
            while (!queryfile.eof()) {
                queryfile >> f1 >> f2;
//                cout<<f1<<" "<<f2<<endl;
                full_list.push_back(make_pair(f1, f2));
            }
            queryfile.close();
            full_list.pop_back();
            cout << "read from queryfile!" << endl;
        }

//        arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);
//        vector<int> arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);

        vector<int> arrival_nodes;

        std::ifstream nodefile;
        string nodefile_name = input_parameters.input_data_dir + "node_"+myversion +configstr+".txt";
        cout<<"node file name:"<<nodefile_name<<endl;
        nodefile.open(nodefile_name, std::ios_base::in);

        if(!nodefile.is_open())
        {
            cout<<"can't load nodefile!"<<endl;
            arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);
            std::ofstream nodesfile_w;
            std::ifstream nodesfile_m;
            nodesfile_w.open(nodefile_name, std::ios_base::out);
            for (int node_i=0;node_i<arrival_nodes.size();node_i++)
            {
                nodesfile_w<<arrival_nodes[node_i]<<endl;
            }
            nodesfile_w.close();
            cout<<"generated nodefile"<<endl;
            nodesfile_m.open(nodefile_name, std::ios_base::in);
            vector<int> arrival_cm;
            while(!nodesfile_m.eof())
            {
                int f3;
                nodesfile_m>>f3;
//            cout<<f3<<endl;
                arrival_cm.push_back(f3);
            }
            nodesfile_m.close();
            cout<<"size:"<<arrival_nodes.size()<<" - "<<arrival_cm.size()<<endl;
            for(int s_i = 0;s_i<arrival_nodes.size();s_i++)
            {
                if(arrival_nodes[s_i]!=arrival_cm[s_i]) cout<<"not equal! "<<s_i<<" "<<arrival_nodes[s_i]<<" - "<<arrival_cm[s_i]<<endl;
            }
            cout<<"Finish Compare!"<<endl<<endl;
        } else{
            while(!nodefile.eof())
            {
                int f3;
                nodefile>>f3;
//            cout<<f3<<endl;
                arrival_nodes.push_back(f3);
            }
            nodefile.close();
            cout<<"read from nodefile!"<<endl;
//            arrival_nodes.erase(arrival_nodes.end()-5,arrival_nodes.end());
            arrival_nodes.pop_back();
            cout<<"full list size:"<<full_list.size()<<" nodes size:"<<arrival_nodes.size()<<endl;
        }

        cout << "full_list made..." << endl;
        cout << "full list size: " << full_list.size() << endl;
        double last_time = 0.0;

        int current_query_threads=0;
        int rand_idx_query=0;
        int rand_idx_update=rand_length;
        int global_start_q_id=0;
        long total_offset=0;
        int total_queries=0;
        int need_opt=0;
        if(multiTestPara.method_name.compare("vtree")==0 && network_name.compare("BJ-old")==0){
            need_opt=1;
        }
        int i;
//        cout<<"init:"<<init_objects<<endl;
        for(i=0;i<init_objects;i++){
//            cout<<i<<endl;

            if(can_estimate)
                gettimeofday(&end, NULL);
            else{
                estimate_mutex.lock();
                gettimeofday(&end, NULL);
                estimate_mutex.unlock();
            }
            long current_time =
                    (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

            pair<double, int> &event = full_list[i];
            int non_object_node = arrival_nodes[i];

            if(need_opt){
                non_object_node%=1270000;
            }

            int start_q_id = global_start_q_id;
            global_start_q_id=(global_start_q_id+1)%num_threads_update;
            rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;

            for(int z = 0; z<num_threads_query;z++) {


                pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                pair<long, pair<int, int> > task = std::make_pair(current_time, node_type_pair);


                // assign insert tasks
                int min_task_size = INT_MAX;

                int pool_index=z * num_threads_update + start_q_id;
                _pool[pool_index]->add_task(task);
                total_object_map[z][non_object_node] = pool_index;
            }
            if(VERIFY){
//                    car_nodes[non_object_node]=1;
                DijkstraKNNInsert(non_object_node, car_nodes);
            }

        }
        if(can_estimate)
            gettimeofday(&global_start, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&global_start, NULL);
            estimate_mutex.unlock();
        }



        struct timeval global_start_2;
        if(can_estimate)
            gettimeofday(&global_start_2, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&global_start_2, NULL);
            estimate_mutex.unlock();
        }

        for (; i < full_list.size(); i++) {
            if(overload_flag) break;
//            cout<<i<<endl;
            if(arrival_nodes[i]==-1) continue;
            pair<double, int> &event = full_list[i];
            long issue_time = floor(event.first * MICROSEC_PER_SEC);
//            cout<<"time: "<<event.first<<" sec"<<endl;
            long current_time=0;
            if(!VERIFY) {

                do {
                    if (issue_time <= current_time) {
//                        if(event.second==QUERY && simulation_time==200) {
//                            if(issue_time + 50 < current_time) {
//                                cout << "current time > issue time + 50" << endl;
//                                cout << "current time: " << current_time << endl;
//                                cout << "issue time : " << issue_time << endl;
//                            }
//                        }
                        break;
                    }
                    if(can_estimate)
                        gettimeofday(&end, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&end, NULL);
                        estimate_mutex.unlock();
                    }
                    current_time =
                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                    if (issue_time > current_time) std::this_thread::sleep_for(std::chrono::microseconds(1));
                } while (true);
//                if(i<init_objects){
//                    issue_time = current_time;
//                }
                issue_time = current_time;
                if(i==init_objects){
                    if(can_estimate)
                        gettimeofday(&global_start, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&global_start, NULL);
                        estimate_mutex.unlock();
                    }


                }
                // start from queue id $start_q_id, we list num_queues_selected consecutive queues to hold random updates
            }

            // if insert
            if (event.second == INSERT) {
                int non_object_node = arrival_nodes[i];

                if(need_opt){
                    non_object_node%=1270000;
                }

//                int start_q_id = global_random_numbers[rand_idx_update] % num_threads_update;
                int start_q_id = global_start_q_id;
                global_start_q_id=(global_start_q_id+1)%num_threads_update;
                rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;
                for(int z = 0; z<num_threads_query;z++) {

                    pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                    pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);

                    // assign insert tasks
                    int min_task_size = INT_MAX;
                    int pool_index=z * num_threads_update + start_q_id;
                    _pool[pool_index]->add_task(task);
                    total_object_map[z][non_object_node] = pool_index;
                }
                if(VERIFY){
//                    car_nodes[non_object_node]=1;
                    DijkstraKNNInsert(non_object_node, car_nodes);
                }
            }
            if (event.second == DELETE) {

                int object_node = arrival_nodes[i];
                if(need_opt){
                    object_node%=1270000;
                }

                for(int z = 0; z<num_threads_query;z++) {
                    int start_q_id = total_object_map[z][object_node];
                    for (int j = start_q_id; j < 1 + start_q_id; j++) {
                        int j_mod = j % num_threads_update;
                        pair<int, int> node_type_pair = std::make_pair(object_node, DELETE);
                        pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
//                        cout<<"delete added to "<<z * num_threads_update + j_mod<<endl;
                        _pool[z * num_threads_update + j_mod]->add_task(task);
                    }
                }
                if(VERIFY){
//                    car_nodes[object_node]=0;
                    DijkstraKNNDelete(object_node, car_nodes);
                }
            }
            if (event.second == QUERY) {

//                cout<<"query "<<endl;
                total_queries++;
                if(can_estimate) {
                    gettimeofday(&end, NULL);
                }
                else{
                    estimate_mutex.lock();
                    gettimeofday(&end, NULL);
                    estimate_mutex.unlock();
                }
                current_time =
                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                total_offset+=current_time-issue_time;

                int query_node=arrival_nodes[i];
                if(need_opt){
                    query_node%=1270000;
                }

                if(VERIFY){
                    vector<KNode> result = DijkstraKNNQuery(k, query_node, mems.dist,
                                                            mems.visited, mems.q, car_nodes);
                    verify_results.push_back(result);

                }
                // put to query tasks
                for (int j = 0; j < num_threads_update; j++) {
                    pair<int, int> node_type_pair = std::make_pair(query_node, QUERY);
                    pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
//                    cout<<"query added to "<<current_query_threads * num_threads_update + j<<endl;
                    if(multiTestPara.method_name.compare("dijk")!=0) {
                    int use_thread_id=-1;
                    long min_cost = INT_MAX;
                    for(int u = current_query_threads;u < current_query_threads+num_threads_query; u++) {
                        int u_mod = u % num_threads_query;
                        long est_cost = _pool[u_mod * num_threads_update + j]->get_est_cost();
                        if(min_cost>est_cost){
                            min_cost = est_cost;
                            use_thread_id = u_mod;
                        }

                    }

                    _pool[use_thread_id * num_threads_update + j]->add_task(task);
                    }
                    else
                        _pool[current_query_threads * num_threads_update + j]->add_task(task);

                }

                current_query_threads=(current_query_threads+1)%num_threads_query;
//                cout<<"query assign cost: "<<clock()-start_1<<endl;
            }
        }
        while(globalThreadVar[threadpool_id][0]->number_of_queries<2){
            std::this_thread::sleep_for(std::chrono::microseconds(1));
        }
        wait_for_finish();
        struct timeval end_2;
        if(can_estimate)
            gettimeofday(&end_2, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&end_2, NULL);
            estimate_mutex.unlock();
        }

        long duration =
                (end_2.tv_sec - global_start_2.tv_sec);
        cout<<"duration: "<<duration<<" secs; fulllist size: "<<full_list.size()<<endl;
        avg_offset = total_offset/total_queries;
        cout<<"avg offset: "<<total_offset/total_queries<<endl;
        if(VERIFY){
            delete[] car_nodes;
            delete_mems(mems);
        }
//        for(int j = 0; j < num_threads_query;j++) {
        for (int i = 0; i < _pool.size(); i++) {
            _pool[i]->set_stop();
            _pool[i]->notify();

        }
        if(!multiTestPara.is_single_aggregate) {
            for (int j = 0; j < num_threads_query; j++) {
                _aggregate_thread[j]->set_stop();
                _aggregate_thread[j]->notify();

            }
        }
        if(multiTestPara.is_single_aggregate){
            _single_aggregate_thread->set_stop();
            _single_aggregate_thread->notify();
        }
        cout << "all set stopped!" << endl;
        Generate_results();
        _needjoin = 1;
        join();
    }

};

#endif //TOAIN_QUERYINGTHREAD_H
