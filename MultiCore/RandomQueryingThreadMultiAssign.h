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
    bool _isfree;
    queue<pair<int, vector<KNode> > > _partial_result_queue;//first query_id; second top-k(k_star) list
    queue<long> issue_time_queue;
    int num_threads_update;
    vector<KNode> knnresult[QUERY_ID_FOLD];
    vector<int> merge_cnt;
    int _stop;
    int k;
    std::mutex _thread_mutex;
    timeval end;
    std::condition_variable cv;



public:

    //构造
    RandomAggregateThread(int copy_id_val, int k_val, int num_threads_update_val) : _isfree(true), _stop(0) {
        copy_id = copy_id_val;
        k = k_val;
        num_threads_update = num_threads_update_val;
        for (int i = 0; i < QUERY_ID_FOLD; i++) merge_cnt.push_back(0);
        _thread = std::thread(&RandomAggregateThread::run, this);
        // _thread.detach(); //放到后台， join是等待线程结束
    }

    void join() {
        _thread.join();
    }

    //是否空闲
    bool isfree() {
        return _isfree;
    }

    void set_stop() {
        _stop = 1;
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
        _thread_mutex.lock();
        issue_time_queue.push(issue_time);
        _partial_result_queue.push(partial_result);
        _thread_mutex.unlock();
        cv.notify_all();
    }

    void notify(){
        cv.notify_all();
    }

    void run() {
        cout << "enter aggregate run!" << endl;
        while (true) {
            if(_partial_result_queue.empty()){
                if(!_stop) {
                    std::unique_lock<std::mutex> aggregate_lock(_thread_mutex);
//                cv.wait(aggregate_lock, []{return true;});
                    cv.wait(aggregate_lock, [this] { return !(_partial_result_queue.empty())||_stop; });
                }


            }

            if (!_partial_result_queue.empty()) {
                _thread_mutex.lock();
                // DO NOT use reference here!!
                pair<int, vector<KNode> > partial_result = _partial_result_queue.front();
                long issue_time = issue_time_queue.front();
                _partial_result_queue.pop();
                issue_time_queue.pop();
                _thread_mutex.unlock();

                int adjust_id = partial_result.first % QUERY_ID_FOLD;
                merge_cnt[adjust_id]++;
                if (merge_cnt[adjust_id] == 1) {
                    knnresult[adjust_id]=partial_result.second;
//                    cout<<"yyyyyyyyyyyyyy: "<<adjust_id<<" "<<knnresult[adjust_id].size()<<endl;
                } else {
//                    cout<<"xxxxxxxxxxxxx!"<<endl;
                    knnresult[adjust_id] = merge_k(knnresult[adjust_id], partial_result.second, k);
                }

                if (merge_cnt[adjust_id] == num_threads_update) {

                    gettimeofday(&end, NULL);
                    long current_time =
                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
                    long response_time = current_time - issue_time;
//                    cout<<"current_time: "<<current_time<<endl;
//                    cout<<"issue_time: "<<issue_time<<endl;

                    globalThreadVar[copy_id]->total_query_time += response_time * 1.0 / MICROSEC_PER_SEC;


                    // cout top-k
                    merge_cnt[adjust_id] = 0;
//                    knnresult[adjust_id].clear();
                    if(VERIFY) {
                        cout<<"start verifying"<<endl;
                        int res = verify_two_results_simple(knnresult[adjust_id], verify_results[adjust_id], k);
                        if(res){
                            cout<<"correct!"<<endl;
                        }
                    }

                    if (is_thresholded) {
                        globalThreadVar[copy_id]->ran_global_locker.lock();
                        globalThreadVar[copy_id]->ran_threshold[adjust_id] = INT_MAX;
                        globalThreadVar[copy_id]->ran_global_locker.unlock();
                    }
                }

            }

//            else if (_stop) {
            _thread_mutex.lock();
            bool stop_cond = _partial_result_queue.empty() && _stop;
            _thread_mutex.unlock();
            if (stop_cond) {
                if (globalThreadVar[copy_id]->stop_run_threads == num_threads_update) {
                    cout << "stopped aggregate threads running" << endl;
                    break;
                }
            }

        }
        cout << endl << endl << "out aggregate" << endl;
    }

};

class RandomThread {
private:
    int copy_id;
    std::thread _thread;
    int _threadid;
    bool _isfree;
    bool _iswait;
    queue<pair<long, pair<int, int> > > _task_queue;
    vector<pair<long, pair<int, int> > > _task_cache_array;
    int _task_index=0;
    int _stop;
    int k;
    mem_struct mems;
    std::mutex _thread_mutex;
    //Algorithm data structure
    int *dijkstra_object_map;
    int query_id = 0;
    RandomAggregateThread *_aggregate_thread;
    std::condition_variable cv;
    timeval end;


public:

    //构造
    RandomThread(int copy_id_val, int id, int k_val, RandomAggregateThread *aggregate_thread_val) : _threadid(id), _isfree(true),
                                                                                                    _stop(0) {
        copy_id = copy_id_val;
        k = k_val;
        _iswait=false;
        allocate_mem(mems, test_n + test_n + 1);
        dijkstra_object_map = new int[test_n + 1];
        memset(dijkstra_object_map, 0, sizeof(int) * (test_n + 1));
        _aggregate_thread = aggregate_thread_val;
        _thread = std::thread(&RandomThread::run, this);
        // _thread.detach(); //放到后台， join是等待线程结束
    }

    void join() {
        _thread.join();
    }

    bool isfree() {
        return _isfree;
    }

    bool iswait(){
        return _iswait;
    }

    void set_stop() {
        _stop = 1;
    }

    void add_task(pair<long, pair<int, int> > &taskids) {
//        if (_isfree) {
        // std::lock_guard<std::mutex>(_locker2);
        _thread_mutex.lock();
        _task_queue.push(taskids);
//            _isfree = false;
        _thread_mutex.unlock();
//        }
        if(_iswait)
            cv.notify_all();
    }

    void add_task_cache(pair<long, pair<int, int> > &taskids) {
        _task_cache_array.push_back(taskids);
        cv.notify_all();

    }

    bool is_task_empty(){
        return _task_index>=_task_cache_array.size();
    }

    void notify(){
        cv.notify_all();
    }
    void run() {
        cout << "enter thread run!" << endl;
        while (true) {

            if(_task_queue.empty()){
                if(!_stop) {

                    std::unique_lock<std::mutex> process_lock(_thread_mutex);
                    _iswait=true;
                    cv.wait(process_lock, [this] { return !(_task_queue.empty())||_stop; });
//                    cv.wait(process_lock, [this] { return true; });
                    _iswait=false;
                }


            }

//            if(is_task_empty()){
//                if(!_stop) {
//
//                    std::unique_lock<std::mutex> process_lock(_thread_mutex);
//                    cv.wait(process_lock, [this] { return !is_task_empty()||_stop; });
////                    cv.wait(process_lock, [this] { return true; });
//                }
//
//
//            }
            if (!_task_queue.empty()) {
//            if(!is_task_empty()){

                _thread_mutex.lock();
                pair<int, int> _taskid = _task_queue.front().second;// extract (node, type) pair
                long _task_time = _task_queue.front().first;
                _task_queue.pop();
                _thread_mutex.unlock();
//                pair<int, int> _taskid = _task_cache_array[_task_index].second;// extract (node, type) pair
//                long _task_time = _task_cache_array[_task_index].first;
//                _task_index++;


                if (_taskid.second == QUERY) {
//                    cout<<"start querying..."<<endl;
                    long start = clock();
                    vector<KNode> kNNs;

                    if (is_thresholded)

                        kNNs = handleQueryThreshold(copy_id, _threadid, multiTestPara.method_name, k, _taskid.first, mems.dist, mems.visited,
                                                    mems.q,
                                                    dijkstra_object_map, globalThreadVar[copy_id]->ran_threshold, query_id);
                    else
                        kNNs = handleQuery(copy_id, _threadid, multiTestPara.method_name, k, _taskid.first, mems.dist, mems.visited, mems.q,
                                           dijkstra_object_map, globalThreadVar[copy_id]->ran_threshold, query_id);

                    if (is_thresholded) {
                        if (kNNs.size() == k) {
                            globalThreadVar[copy_id]->ran_global_locker.lock();

                            if (kNNs[k - 1].dis < globalThreadVar[copy_id]->ran_threshold[query_id])
                                globalThreadVar[copy_id]->ran_threshold[query_id] = kNNs[k - 1].dis;
                            globalThreadVar[copy_id]->ran_global_locker.unlock();
                        }
                    }
                    if(num_threads_update>1) {
                        pair<int, vector<KNode> > partial_res = make_pair(query_id, kNNs);
                        query_id = (query_id + 1) % QUERY_ID_FOLD;
                        _aggregate_thread->add_task(_task_time, partial_res);
                    }
                    else{
                        gettimeofday(&end, NULL);
                        long current_time =
                                (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
                        long response_time = current_time - _task_time;
//
                        globalThreadVar[copy_id]->total_query_time += response_time * 1.0 / MICROSEC_PER_SEC;

                    }
//                    cout<<"aggregate task added!"<<endl;

                    // put partial querying result to somewhere
                }
                if (_taskid.second == INSERT) {
                    handleInsert(copy_id, _threadid, multiTestPara.method_name, k, _taskid.first, dijkstra_object_map, mems);

                }
                if (_taskid.second == DELETE) {
                    handleDelete(copy_id, _threadid, multiTestPara.method_name, k, _taskid.first, dijkstra_object_map, mems);

                }

            }

            _thread_mutex.lock();
            bool cond_stop = _task_queue.empty() && _stop;
            _thread_mutex.unlock();
//            _thread_mutex.lock();
//            bool cond_stop = is_task_empty() && _stop;
//            _thread_mutex.unlock();
            if (cond_stop) {
//                cout << "stopped running" << endl;
                globalThreadVar[copy_id]->ran_global_locker.lock();
                globalThreadVar[copy_id]->stop_run_threads++;
                globalThreadVar[copy_id]->ran_global_locker.unlock();
                gettimeofday(&end, NULL);
                long current_time =
                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
                cout<<"copyid: "<<copy_id<<" threadid: "<<_threadid<<" ending time: "<<current_time<<endl;
                break;
            }

        }
//        cout << endl << endl << "out 2" << endl;
    }
};



class RandomThreadPool {
private:
    int threadpool_id;
    int begin_node;
    int end_node;
    vector<RandomThread *> _pool;
     std::mutex pool_mutex;
    std::thread* _main_thread;
    RandomAggregateThread **_aggregate_thread;

    RandomAggregateThread * _single_aggregate_thread;
    int num_threads_update;
    int _needjoin;
    double alpha;
    int test_n;
    double fail_p;
    int k;
    // stores the part where the object locates
    vector<int> total_object_map;
    vector<int> distribute_sizes;

    double query_rate;
    double insert_rate;
    double delete_rate;
    int num_threads_query;
    int simulation_time;
    vector<std::pair<double, int> > full_list;
    vector<int> arrival_nodes;

    queue<std::pair<double, int> > arrival_tasks_queue;
    queue<int> arrival_nodes_queue;
    int num_dispatch_thread;
    int global_init_objects;


    // Algorithm data structure
//    vector<int *> dijkstra_object_map_vec;
public:
    RandomThreadPool(int threadpool_id_val, int begin_node_val, int end_node_val, int num_threads_query_val, int num_threads_update_val, double alpha_val, int k_val, double fail_p_val,
                     int test_n_val, double query_rate_val, double insert_rate_val, double delete_rate_val, int simulation_time_val, int num_dispatch_thread_val) {
        cout << "constructing RandomThreadPool..." << endl;
        num_dispatch_thread=num_dispatch_thread_val;
        threadpool_id = threadpool_id_val;
        begin_node = begin_node_val;
        end_node = end_node_val;
        simulation_time = simulation_time_val;
        num_threads_query=num_threads_query_val;
        query_rate = query_rate_val;
        insert_rate = insert_rate_val;
        delete_rate = delete_rate_val;
        test_n = test_n_val;
        num_threads_update = num_threads_update_val;
        alpha = alpha_val;
        fail_p = fail_p_val;
        k = k_val;

        globalThreadVar = new GlobalThreadVar*[num_threads_query];
        if(!is_single_aggregate)
            _aggregate_thread = new RandomAggregateThread*[num_threads_query];
        if(is_single_aggregate){
            _single_aggregate_thread=new RandomAggregateThread(0, k, num_threads_update);
        }
        int k_star = compute_k_star(k, num_threads_update, alpha, fail_p);
        for(int j =0;j<num_threads_query;j++){
            globalThreadVar[j]=new GlobalThreadVar();
            globalThreadVar[j]->ran_threshold.clear();
            for (int i = 0; i <= QUERY_ID_FOLD; i++) globalThreadVar[j]->ran_threshold.push_back(INT_MAX);
//            int i = num_threads_update;

            cout << "k_star: " << k_star << endl;
            if(!is_single_aggregate)
                _aggregate_thread[j] = new RandomAggregateThread(j, k, num_threads_update);
            for(int i=0;i<num_threads_update;i++) {
                if(!is_single_aggregate) {
                    RandomThread *t = new RandomThread(j, i, k_star, _aggregate_thread[j]);
                    _pool.push_back(t);
                }
                else{
                    RandomThread *t = new RandomThread(j, i, k_star, _single_aggregate_thread);
                    _pool.push_back(t);
                }
            }
        }

        _needjoin = 0;
        init();

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
        mem_struct mems;
        int* car_nodes;
        if(VERIFY){
            allocate_mem(mems, test_n+1);
            car_nodes = new int[test_n+1];
            memset(car_nodes, 0, sizeof(int)*(test_n+1));
        }
        global_init_objects=1000;
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
        for (int i = 0; i <= test_n; i++) {
            total_object_map.push_back(0);
            distribute_sizes.push_back(0);
        }
        for (int i = 0; i < test_n; i++) {
            if (!total_object_map[i]) {
                non_object_list.push_back(i);
            }
        }
        int init_objects = 1000;

        for (int i = 0; i < init_objects; i++) {
            full_list.push_back(make_pair(0.0, INSERT));

        }
        vector<std::pair<double, int> > append_list = make_online_query_update_list(query_rate, insert_rate,
                                                                                    delete_rate,
                                                                                    simulation_time);
        for (pair<double, int> &item : append_list)
            full_list.push_back(item);

        arrival_nodes = generate_arrival_nodes(full_list, begin_node, end_node);
        cout << "full_list made..." << endl;
        cout << "full list size: " << full_list.size() << endl;

        for(std::pair<double, int> task: full_list){
            arrival_tasks_queue.push(task);
        }
        for(int node: arrival_nodes){
            arrival_nodes_queue.push(node);
        }

        _main_thread = new std::thread[num_dispatch_thread];
        for(int j=0; j<num_dispatch_thread; j++) {

            _main_thread[j] = std::thread(&RandomThreadPool::run, this, j);
        }
        for(int j=0; j<num_dispatch_thread; j++) {

            _main_thread[j].join();
        }
        delete[] _main_thread;

        if(VERIFY){
            delete[] car_nodes;
            delete_mems(mems);
        }
//        for(int j = 0; j < num_threads_query;j++) {
        for (int i = 0; i < _pool.size(); i++) {
            _pool[i]->set_stop();
            _pool[i]->notify();
        }
        if(!is_single_aggregate) {
            for (int j = 0; j < num_threads_query; j++) {
                _aggregate_thread[j]->set_stop();
                _aggregate_thread[j]->notify();

            }
        }
        if(is_single_aggregate){
            _single_aggregate_thread->set_stop();
        }
        _needjoin = 1;
        cout << "all set stopped!" << endl;

    }

    void join() {
        cout << "start joining" << endl;
//        _main_thread.join();
        for (int i = 0; i < _pool.size(); i++) {
            _pool[i]->join();
        }
        if(!is_single_aggregate) {
            for (int j = 0; j < num_threads_query; j++) {
                _aggregate_thread[j]->join();
            }
        }
        if(is_single_aggregate){
            _single_aggregate_thread->join();
        }
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

//        _main_thread = std::thread(&RandomThreadPool::run, this);
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
                int index = global_random_numbers[rand_idx_update] % non_object_list.size();
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
                int index = global_random_numbers[rand_idx_update] % object_list.size();
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
                int query_node;
                do {
                    query_node = global_random_numbers[rand_idx_query];
                    rand_idx_query=(rand_idx_query+1)%rand_length;
                    query_node = query_node % (test_n-1)+1;

                } while (vStart[query_node] == 0);
                arrival_node_list.push_back(query_node);

            }

        }
        return arrival_node_list;

    }

    void run(int start_id) {

        struct timeval start;
        struct timeval end;
        int init_objects = global_init_objects/num_dispatch_thread;
        double last_time = 0.0;
        gettimeofday(&global_start, NULL);
        struct timeval global_start_2;
        gettimeofday(&global_start_2, NULL);
        int current_query_threads=0;
        int rand_idx_query=0;
        int rand_idx_update=rand_length;
        int global_start_q_id=0;
        long total_offset=0;
        int total_queries=0;
        long issue_time=0;
        long pre_issue_time = 0;
        long current_time=0;
        while(true) {
            if(arrival_nodes.empty()) break;
            pre_issue_time=issue_time;
            pool_mutex.lock();
            int arrival_node = arrival_nodes_queue.front();
            arrival_nodes_queue.pop();
            std::pair<double, int> event = arrival_tasks_queue.front();
            arrival_tasks_queue.pop();
            pool_mutex.unlock();

            if(arrival_node==-1) continue;
                issue_time = floor(event.first * MICROSEC_PER_SEC);
//            cout<<"time: "<<event.first<<" sec"<<endl;

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
                    gettimeofday(&end, NULL);
                    current_time =
                            (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                    if (issue_time > current_time) std::this_thread::sleep_for(std::chrono::microseconds(1));
                } while (true);
//                if(i<init_objects){
//                    issue_time = current_time;
//                }
                if(pre_issue_time==0 && issue_time > 0){
                    gettimeofday(&global_start, NULL);

                }
                // start from queue id $start_q_id, we list num_queues_selected consecutive queues to hold random updates
            }


            // if insert
            if (event.second == INSERT) {

                int non_object_node = arrival_node;


//                int start_q_id = global_random_numbers[rand_idx_update] % num_threads_update;
                int start_q_id = global_start_q_id;
                global_start_q_id=(global_start_q_id+1)%num_threads_update;
                rand_idx_update=(rand_idx_update+1)%rand_length+rand_length;
//                object_list.push_back(non_object_node);

                total_object_map[non_object_node] = start_q_id;
//                cout<<"here"<<endl;
//                distribute_sizes[non_object_node] = num_queues_selected;
                pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                pair<long, pair<int, int> > task = std::make_pair(current_time, node_type_pair);
                for (int j = start_q_id; j < 1 + start_q_id; j++) {
                    int j_mod = j % num_threads_update;

                    // assign insert tasks

                    for(int z = 0; z<num_threads_query;z++) {
//                        cout<<num_threads_update<<" "<<z<<" "<<_pool.size()<<" "<<z * num_threads_update + j_mod<<endl;
                        _pool[z * num_threads_update + j_mod]->add_task(task);
                    }
                }
//                if(VERIFY){
////                    car_nodes[non_object_node]=1;
//                    DijkstraKNNInsert(non_object_node, car_nodes);
//                }
//                cout<<"insert assign cost : "<<clock()-start_1<<endl;
//                cout<<"end insert"<<endl;
            }
            if (event.second == DELETE) {
                int object_node = arrival_node;
                int start_q_id = total_object_map[object_node];
                pair<int, int> node_type_pair = std::make_pair(object_node, DELETE);
                pair<long, pair<int, int> > task = std::make_pair(current_time, node_type_pair);
                for (int j = start_q_id; j < distribute_sizes[object_node] + start_q_id; j++) {
                    int j_mod = j % num_threads_update;

                    for(int z = 0; z<num_threads_query;z++)
                        _pool[z*num_threads_update+j_mod]->add_task(task);
                }
//                if(VERIFY){
////                    car_nodes[object_node]=0;
//                    DijkstraKNNDelete(object_node, car_nodes);
//                }

//                cout<<"delete assign cost: "<<clock()-start_1<<endl;

            }
            if (event.second == QUERY) {
                total_queries++;

                gettimeofday(&end, NULL);
                current_time =
                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                total_offset+=current_time-issue_time;


                globalThreadVar[current_query_threads]->number_of_queries++;

//                do {
//                    query_node = global_random_numbers[rand_idx_query];
//                    rand_idx_query=(rand_idx_query+1)%rand_length;
//                    query_node = query_node % (test_n-1)+1;
//
//                } while (vStart[query_node] == 0);
                int query_node=arrival_node;

//                if(VERIFY){
//                    vector<KNode> result = DijkstraKNNQuery(k, query_node, mems.dist,
//                                                            mems.visited, mems.q, car_nodes);
//                    verify_results.push_back(result);
//
//                }
                // put to query tasks
//                gettimeofday(&end, NULL);
//                current_time =
//                        (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;

                for (int j = 0; j < num_threads_update; j++) {
                    pair<int, int> node_type_pair = std::make_pair(query_node, QUERY);
                    pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
                    _pool[current_query_threads * num_threads_update + j]->add_task(task);


                }

                current_query_threads=(current_query_threads+1)%num_threads_query;
//                cout<<"query assign cost: "<<clock()-start_1<<endl;
            }

        }
        struct timeval end_2;
        gettimeofday(&end_2, NULL);
        long duration =
                (end_2.tv_sec - global_start_2.tv_sec);
        cout<<"duration: "<<duration<<" secs; fulllist size: "<<full_list.size()<<endl;
        cout<<"avg offset: "<<total_offset/total_queries<<endl;

    }

};

#endif //TOAIN_QUERYINGTHREAD_H
