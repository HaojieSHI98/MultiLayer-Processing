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



public:

    //构造
    RandomAggregateThread(int copy_id_val, int k_val, int num_threads_update_val) : is_free(true), thread_stop(0) {
        copy_id = copy_id_val;
        k = k_val;
        num_threads_update = num_threads_update_val;
        for (int i = 0; i < QUERY_ID_FOLD; i++) merge_cnt.push_back(0);
        _thread = std::thread(&RandomAggregateThread::run, this);
        // working_thread.detach(); //放到后台， join是等待线程结束
    }

    void join() {
        if(_thread.joinable())
            _thread.join();
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
        thread_mutex.unlock();
//        if(is_wait)
            cv.notify_all();
//        cout<<"aggregate task size: "<<partial_result_queue.size()<<endl;
    }

    void notify(){
        cv.notify_all();
    }

    void run() {
        cout << "enter aggregate run!" << endl;
        while (true) {
            if(overload_flag) break;
            if(is_simulation && thread_stop){
                break;
            }
//            if(!thread_stop)
//             {
                 if(!is_wait) {
                     if (partial_result_queue.empty()) {
                         std::unique_lock<std::mutex> aggregate_lock(thread_mutex);
//                cv.wait(aggregate_lock, []{return true;});
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
                thread_mutex.unlock();

                int adjust_id = partial_result.first % QUERY_ID_FOLD;
                merge_cnt[adjust_id]++;
                if (merge_cnt[adjust_id] == 1) {
                    knn_result[adjust_id]=partial_result.second;
//                    cout<<"yyyyyyyyyyyyyy: "<<adjust_id<<" "<<knn_result[adjust_id].size()<<endl;
                } else {
//                    cout<<"xxxxxxxxxxxxx!"<<endl;
                    knn_result[adjust_id] = merge_k(knn_result[adjust_id], partial_result.second, k);
                }

                if (merge_cnt[adjust_id] == num_threads_update) {

                    if(can_estimate)
                        gettimeofday(&end, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&end, NULL);
                        estimate_mutex.unlock();
                    }
                    long current_time =
                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
                    long response_time = current_time - issue_time;
                    response_time_list.push_back(response_time);
//                    if(response_time>100*MICROSEC_PER_SEC){
//                        cout<<"overloaded!"<<endl;
//                        overload_flag = 1;
//                        break;
//                    }
//                    cout<<"current_time: "<<current_time<<endl;
//                    cout<<"issue_time: "<<issue_time<<endl;

                    globalThreadVar[copy_id]->total_query_time += response_time * 1.0 / MICROSEC_PER_SEC;
                    globalThreadVar[copy_id]->number_of_queries++;


                    // cout top-k
                    merge_cnt[adjust_id] = 0;
                    if(VERIFY) {
                        cout<<"start verifying"<<endl;
                        int res = verify_two_results_simple(knn_result[adjust_id], verify_results[adjust_id], k);
                        if(res){
                            cout<<"correct!"<<endl;
                        }
                    }

                    if (multiTestPara.is_thresholded) {
                        globalThreadVar[copy_id]->ran_global_locker.lock();
                        globalThreadVar[copy_id]->ran_threshold[adjust_id] = INT_MAX;
                        globalThreadVar[copy_id]->ran_global_locker.unlock();
                    }
                }

            }

//            else if (thread_stop) {
            thread_mutex.lock();
            bool stop_cond = partial_result_queue.empty() && thread_stop;
            thread_mutex.unlock();
            if (stop_cond) {
                if(!multiTestPara.is_single_aggregate) {
                    if (globalThreadVar[copy_id]->stop_run_threads == num_threads_update) {
                        cout << "stopped aggregate threads running" << endl;
                        break;
                    }
                }
                else{
                    if (globalThreadVar[copy_id]->stop_run_threads == num_threads_update*multiTestPara.num_threads_query) {
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


public:

    //构造
    RandomThread(int copy_id_val, int id, int k_val,int query_cost,int insert_cost,int delete_cost, RandomAggregateThread *aggregate_thread_val) : thread_id(id), is_free(true),
                                                                                                    thread_stop(0) {
        copy_id = copy_id_val;
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
        cout<<"query_lost:"<<last_query_cost<<"insert_lost:"<<last_insert_cost<<"delete_lost:"<<last_delete_cost<<endl;

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
        cout << "enter thread run!" << endl;
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
                                                    dijkstra_object_map, globalThreadVar[copy_id]->ran_threshold,
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
                                           dijkstra_object_map, globalThreadVar[copy_id]->ran_threshold, query_id,
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

                    }
                    if (multiTestPara.is_thresholded) {
                        if (kNNs.size() == k) {
                            globalThreadVar[copy_id]->ran_global_locker.lock();

                            if (kNNs[k - 1].dis < globalThreadVar[copy_id]->ran_threshold[query_id])
                                globalThreadVar[copy_id]->ran_threshold[query_id] = kNNs[k - 1].dis;
                            globalThreadVar[copy_id]->ran_global_locker.unlock();
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
                    globalThreadVar[copy_id]->ran_global_locker.lock();
                    globalThreadVar[copy_id]->stop_run_threads++;
                    globalThreadVar[copy_id]->ran_global_locker.unlock();
                }
                else {
                    globalThreadVar[0]->ran_global_locker.lock();
                    globalThreadVar[0]->stop_run_threads++;
                    globalThreadVar[0]->ran_global_locker.unlock();
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


    // Algorithm data structure
//    vector<int *> dijkstra_object_map_vec;
public:
    RandomThreadPool(int threadpool_id_val, int begin_node_val, int end_node_val, int num_threads_query_val, int num_threads_update_val, double alpha_val, int k_val, double fail_p_val,
                     int test_n_val, double query_rate_val, double insert_rate_val, double delete_rate_val, int simulation_time_val,int query_cost,int insert_cost,int delete_cost) {
        cout << "constructing RandomThreadPool..." << endl;
        if(overload_flag)
            overload_flag=0;
        total_queries_plan=0;
        total_updates_plan=0;
        total_queries_finished=0;
        total_updates_finished=0;
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
        alpha = alpha_val;
        fail_p = fail_p_val;
        k = k_val;
        init();
        cout<<"num_threads_query: "<<num_threads_query<<endl;
        globalThreadVar = new GlobalThreadVar*[num_threads_query];
        int k_star = compute_k_star(k, num_threads_update, alpha, fail_p);
        for(int j =0;j<num_threads_query;j++) {
            globalThreadVar[j] = new GlobalThreadVar();
            globalThreadVar[j]->ran_threshold.clear();
            for (int i = 0; i <= QUERY_ID_FOLD; i++) globalThreadVar[j]->ran_threshold.push_back(INT_MAX);
        }
        if(!multiTestPara.is_single_aggregate)
            _aggregate_thread = new RandomAggregateThread*[num_threads_query];
        if(multiTestPara.is_single_aggregate){
            _single_aggregate_thread=new RandomAggregateThread(0, k, num_threads_update);
        }
        cout << "k_star: " << k_star << endl;
        for(int j =0;j<num_threads_query;j++){

            if(!multiTestPara.is_single_aggregate)
                _aggregate_thread[j] = new RandomAggregateThread(j, k, num_threads_update);
            for(int i=0;i<num_threads_update;i++) {
                if(!multiTestPara.is_single_aggregate) {
                    RandomThread *t = new RandomThread(j, i, k_star,query_cost,insert_cost,delete_cost, _aggregate_thread[j]);
                    _pool.push_back(t);
                }
                else{
                    RandomThread *t = new RandomThread(j, i, k_star,query_cost,insert_cost,delete_cost, _single_aggregate_thread);
                    _pool.push_back(t);
                }
            }
        }

        _needjoin = 0;
    }

    //释放线程池
    ~RandomThreadPool() {
        cout << "hi, RandomThreadPooladPool()" << endl;
        // releasePool();
        // for(int i = 0;i < _pool.size(); ++i){
        //     delete _pool[i];
        // }

    }

    void init() {
        cout << "start init() function ..." << endl;

    }

    void join() {
        cout << "start joining" << endl;

        if(!multiTestPara.is_single_aggregate) {
            for (int j = 0; j < num_threads_query; j++) {
                _aggregate_thread[j]->join();
            }
        }
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
        for (int i = 0; i < multiTestPara.num_threads_query; i++) {
            total_queries_finished += globalThreadVar[i]->number_of_queries;
        }
        query_finish_rate = total_queries_finished * 1.0 / total_queries_plan;

        cout << "finish joining threads" << endl;

        if(_main_thread.joinable())
            _main_thread.join();

        cout << "finish joining main thread" << endl;

    }

    int isNeedJoin() {
        return _needjoin;
    }

    void releasePool() {
//        for(int j=0;j<num_threads_query;j++) {
        for (int i = 0; i < _pool.size(); ++i) {
            delete _pool[i];
//            }
        }
    }

    void start() {

        _main_thread = std::thread(&RandomThreadPool::run, this);
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

        std::ifstream queryfile;
        queryfile.open(input_parameters.input_data_dir + "query_" +std::to_string(query_rate)+"_"+std::to_string(insert_rate)+"_"+
                               std::to_string(delete_rate)+"_"+std::to_string(simulation_time)+".txt", std::ios_base::in);
        double f1;
        int f2;
        if(!queryfile.is_open())
        {
            cout<<"can't load queryfile!"<<endl;
            vector<std::pair<double, int> > append_list = make_online_query_update_list(query_rate, insert_rate,
                                                                                         delete_rate,
                                                                                         simulation_time);
            std::ofstream queryfile_w;
            queryfile_w.open(input_parameters.input_data_dir + "query_" +std::to_string(query_rate)+"_"+std::to_string(insert_rate)+"_"+
                           std::to_string(delete_rate)+"_"+std::to_string(simulation_time)+".txt", std::ios_base::out);
            for (pair<double, int> &item : append_list)
             {
                 full_list.push_back(item);
                 queryfile_w<<item.first<<" "<<item.second<<endl;
             }
            queryfile_w.close();
            cout<<"write to queryfile!"<<endl;
        }
        else{
            while(!queryfile.eof())
            {
                queryfile>>f1>>f2;
//                cout<<f1<<" "<<f2<<endl;
                full_list.push_back(make_pair(f1,f2));
            }
            queryfile.close();
            cout<<"read from queryfile!"<<endl;
        }

        vector<int> arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);

//        vector<int> arrival_nodes;
//        std::ifstream nodefile;
//        nodefile.open(input_parameters.input_data_dir + "node_" +std::to_string(query_rate)+"_"+std::to_string(insert_rate)+"_"+
//                      std::to_string(delete_rate)+"_"+std::to_string(simulation_time)+".txt", std::ios_base::in);
//        int f3;
//        if(!nodefile.is_open())
//        {
//            cout<<"can't load nodefile!"<<endl;
//            arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);
//            std::ofstream nodesfile_w;
//            nodesfile_w.open(input_parameters.input_data_dir + "node_" +std::to_string(query_rate)+"_"+std::to_string(insert_rate)+"_"+
//                             std::to_string(delete_rate)+"_"+std::to_string(simulation_time)+".txt", std::ios_base::out);
//             for (int node_i=0;node_i<arrival_nodes.size();node_i++)
//             {
//                 nodesfile_w<<arrival_nodes[node_i]<<endl;
//             }
//            nodesfile_w.close();
//        } else{
//            while(!nodefile.eof())
//            {
//                nodefile>>f3;
////            cout<<f3<<endl;
//                arrival_nodes.push_back(f3);
//            }
//            nodefile.close();
//            cout<<"read from nodefile!"<<endl;
//        }

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
                if(i<init_objects){
                    issue_time = current_time;
                }
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
//            cout<<"get here"<<endl;
//            if(!VERIFY) {
//                long current_time;
////                do {
//                    gettimeofday(&end, NULL);
//                    current_time =
//                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
//                    if (issue_time <= current_time) {
//                        if(event.second==QUERY) {
//                            cout << "current time >= issue time" << endl;
//                            cout << "current time: " << current_time << endl;
//                            cout << "issue time : " << issue_time << endl;
//                        }
//                    }
//                    else std::this_thread::sleep_for(std::chrono::microseconds(issue_time-current_time));
////                } while (true);
//                // start from queue id $start_q_id, we list num_queues_selected consecutive queues to hold random updates
//                if(i<init_objects){
//                    issue_time = current_time;
//                }
//                if(i==init_objects){
//                    gettimeofday(&global_start, NULL);
//
//                }
//            }


            // if insert
            if (event.second == INSERT) {




//                long start_1 = clock();
//                int index = global_random_numbers[rand_idx_update] % non_object_list.size();
//                rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;
//                int non_object_node = non_object_list[index];
//                int last_index = non_object_list.size() - 1;
//
//                int tmp = non_object_list[index];
//                non_object_list[index] = non_object_list[last_index];
//                non_object_list[last_index] = tmp;
//                non_object_list.pop_back();
//
//                object_list.push_back(non_object_node);
//                cout<<"INSERT "<< endl;
                int non_object_node = arrival_nodes[i];

                if(need_opt){
                    non_object_node%=1270000;
                }



//                int start_q_id = global_random_numbers[rand_idx_update] % num_threads_update;
                int start_q_id = global_start_q_id;
                global_start_q_id=(global_start_q_id+1)%num_threads_update;
                rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;
//                object_list.push_back(non_object_node);


//                cout<<"here"<<endl;
//                distribute_sizes[non_object_node] = num_queues_selected;

                for(int z = 0; z<num_threads_query;z++) {


                    pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                    pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);


                    // assign insert tasks
                    int min_task_size = INT_MAX;
//                    int pool_index = -1;
//                    for (int j = start_q_id; j < start_q_id+num_threads_update; j++) {
//                        int j_mod = j % num_threads_update;
//
//                        if(_pool[z * num_threads_update + j_mod]->get_task_size()<min_task_size){
//                            pool_index=z * num_threads_update + j_mod;
//                            min_task_size=_pool[z * num_threads_update + j_mod]->get_task_size();
//
//                        }
////                        cout<<"insert added to : "<<z * num_threads_update + j_mod<<endl;
////                        cout<<num_threads_update<<" "<<z<<" "<<_pool.size()<<" "<<z * num_threads_update + j_mod<<endl;
//
//                    }
                    int pool_index=z * num_threads_update + start_q_id;
                    _pool[pool_index]->add_task(task);
                    total_object_map[z][non_object_node] = pool_index;
                }
                if(VERIFY){
//                    car_nodes[non_object_node]=1;
                    DijkstraKNNInsert(non_object_node, car_nodes);
                }
//                cout<<"insert assign cost : "<<clock()-start_1<<endl;
//                cout<<"end insert"<<endl;
            }
            if (event.second == DELETE) {


//                long start_1 = clock();
//                int index = global_random_numbers[rand_idx_update] % object_list.size();
//                rand_idx_update=(rand_idx_update+1)%rand_length;
//                int object_node = object_list[index];
//                int last_index = object_list.size() - 1;
//                int tmp = object_list[index];
//                object_list[index] = object_list[last_index];
//                object_list[last_index] = tmp;
//                object_list.pop_back();
//
//                non_object_list.push_back(object_node);
//                cout<<"DELETE "<<endl;
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

//                cout<<"delete assign cost: "<<clock()-start_1<<endl;

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
//                gettimeofday(&end, NULL);
//                current_time =
//                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

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
        while(globalThreadVar[0]->number_of_queries<2){
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
        _needjoin = 1;
        cout << "all set stopped!" << endl;
    }

};

#endif //TOAIN_QUERYINGTHREAD_H
