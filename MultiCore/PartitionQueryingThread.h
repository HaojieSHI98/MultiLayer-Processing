//
// Created by lostrong on 18/3/19.
//

#ifndef TOAIN_PARTITIONQUERYINGTHREAD_H
#define TOAIN_PARTITIONQUERYINGTHREAD_H


//linux, g++ -std=c++14 -o t *.cpp -pthread
#include <queue>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>
//#include "../multicore/DijkstraKNN.h"
#include "FunctionTemplate.h"
#include "../util/MemoryAllocation.h"
#include "../graph/CellNode.h"
#include "../queue/GeneralQueueSimulator.h"
#include "../MultiCore/Simulation.h"
#include "../MultiCore/ComputeKStar.h"
#include "MultiCoreVariables.h"
#include "util.h"
#include <sys/time.h>

std::mutex global_locker;
std::mutex global_locker_2;
int partition_stop_run_threads = 0;
vector<int> threshold;

class PartitionAggregateThread {
private:
    std::thread _thread;
    bool _isfree;
    queue<pair<int, vector<KNode> > > _partial_result_queue;//first query_id; second top-k(k_star) list
    queue<long> issue_time_queue;
    int num_threads_update;
    vector<KNode> knnresult[QUERY_ID_FOLD];
    vector<int> merge_cnt;
    int _stop;
    int k;
    std::mutex _thread_locker;


public:

    //构造
    PartitionAggregateThread(int k_val, int num_threads_update_val) : _isfree(true), _stop(0) {
        k = k_val;
        num_threads_update = num_threads_update_val;
        for (int i = 0; i < QUERY_ID_FOLD; i++) {
            merge_cnt.push_back(0);

        }
        _thread = std::thread(&PartitionAggregateThread::run, this);
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
//        cout<<"size 1: "<<S1.size()<<" "<<"size 2: "<<S2.size()<<endl;
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

    //添加任务
    void add_task(long issue_time, pair<int, vector<KNode> >& partial_result) {
//        if (_isfree) {
        // std::lock_guard<std::mutex>(_locker2);
//        cout<<"here 1"<<endl;
        _thread_locker.lock();
//        cout<<"here2"<<endl;
        issue_time_queue.push(issue_time);
        _partial_result_queue.push(partial_result);
//            _isfree = false;
        _thread_locker.unlock();
//        }
    }

    //如果有任务则执行任务，否则自旋
    void run() {
        cout << "enter aggregate run!" << endl;
        while (true) {
//            cout<<"hi"<<endl;
//            std::this_thread::sleep_for(std::chrono::microseconds(1));

            if (!_partial_result_queue.empty()) {

//                    cout<<"enter aggregate thread real logic"<<endl;

//                cout<<"size: "<<_partial_result_queue.size()<<endl;

                _thread_locker.lock();
                pair<int, vector<KNode> > partial_result = _partial_result_queue.front();
                long issue_time = issue_time_queue.front();

                _partial_result_queue.pop();
                issue_time_queue.pop();
//                cout<<"partial result size: "<<partial_result.second.size()<<endl;
//                _isfree = true;
                _thread_locker.unlock();
//                _isfree = false;
                //do something

                int adjust_id = partial_result.first % QUERY_ID_FOLD;
                merge_cnt[adjust_id]++;
//                cout<<"begin merging..."<<endl;
//                cout<<"size:"<<partial_result.second.size()<<endl;
                if(merge_cnt[adjust_id]==1) {
                    knnresult[adjust_id] = partial_result.second;
//                    cout<<"yyyyyyyyyyyyyy: "<<adjust_id<<" "<<knnresult[adjust_id].size()<<endl;
                }
                else {
//                    cout<<"xxxxxxxxxxxxx!"<<endl;
                    knnresult[adjust_id] = merge_k(knnresult[adjust_id], partial_result.second, k);
                }
//                cout<<"after merging..."<<endl;
                timeval end;
                if (merge_cnt[adjust_id] == num_threads_update) {
                    // cout top-k
                    merge_cnt[adjust_id] = 0;
                    gettimeofday(&end, NULL);
                    long current_time = (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
                    long response_time = current_time-issue_time;
//                    cout<<"current_time: "<<current_time<<endl;
//                    cout<<"issue_time: "<<issue_time<<endl;
                    total_query_time+=response_time*1.0/MICROSEC_PER_SEC;


                    if(VERIFY) {
                        cout<<"start verifying"<<endl;
                        int res = verify_two_results_simple(knnresult[adjust_id], verify_results[adjust_id], k);
                        if(res){
                            cout<<"correct!"<<endl;
                        }
                    }
//                    knnresult[adjust_id].clear();
                    if(is_thresholded) {
                        global_locker.lock();
                        threshold[adjust_id] = INT_MAX;
                        global_locker.unlock();
                    }
                }

            }

//            else if (_stop) {
            _thread_locker.lock();
            bool stop_cond=_partial_result_queue.empty() && _stop;
            _thread_locker.unlock();
            if(stop_cond){
//                cout<<"num_threads_update: "<<num_threads_update<<endl;
//                cout<<"partition_stop_run_threads: "<<partition_stop_run_threads<<endl;
                if (partition_stop_run_threads == num_threads_update) {
                    cout << "stopped aggregate threads running" << endl;
                    break;
                }
            }

        }
//        cout << endl << endl << "out aggregate" << endl;
    }

};

class PartitionThread {
private:
    std::thread _thread;
    int _threadid;
    bool _isfree;
    queue<pair<long, pair<int, int> > > _task_queue;
    int _stop;
    int k;
    mem_struct mems;
    std::mutex _thread_locker;
    //Algorithm data structure
    int *dijkstra_object_map;
    int query_id = 0;
    vector<int> part_maps;
    PartitionAggregateThread *_aggregate_thread;


public:

    //构造
    PartitionThread(int id, vector<int> &part_maps_val, int k_val, PartitionAggregateThread *aggregate_thread_val)
            : _threadid(id), _isfree(true), _stop(0) {
        k = k_val;
        part_maps = part_maps_val;
        allocate_mem(mems, test_n + test_n + 1);
        dijkstra_object_map = new int[test_n + 1];
        memset(dijkstra_object_map, 0, sizeof(int) * (test_n + 1));
//        if(method_name.compare("toain")==0){
//
//                for(int i = 0;i<test_n;i++){
//                    hier_local_knn_arr_multi[id][i].clear();
//                }
//
//        }
//        if(method_name.compare("vtree")==0){
//
//        }

        _aggregate_thread = aggregate_thread_val;
        _thread = std::thread(&PartitionThread::run, this);
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

    //添加任务
    void add_task(pair<long, pair<int, int> >& taskid) {
//        if (_isfree) {
        // std::lock_guard<std::mutex>(_locker2);
        _thread_locker.lock();
        _task_queue.push(taskid);
//        cout<<"_task_queue size: "<<_task_queue.size()<<endl;
//            _isfree = false;
        _thread_locker.unlock();
//        cout<<"task size: "<<_task_queue.size()<<endl;
//        }
    }

    //如果有任务则执行任务，否则自旋
    void run() {
        cout << "enter thread run!" << endl;
        int first=1;
        while (true) {
//            std::this_thread::sleep_for(std::chrono::nanoseconds(1));

//            if(first)
//            {
//                std::this_thread::sleep_for(std::chrono::microseconds(1));
//                first=0;
//            }
//            if (!_task_queue.empty()) {
//            cout<<"_task_queue.size(): "<<_task_queue.size()<<endl;
            if (!_task_queue.empty()) {

//                if(_threadid==0){
//                    cout<<"enter partitionthread real logic"<<endl;
//                }
//                cout<<_threadid<<endl;
//                cout<<"xxxxxxxxxxxxx!!!"<<endl;
                _thread_locker.lock();
//                _isfree = false;

                //do something
                pair<int, int> _taskid = _task_queue.front().second;// extract (node, type) pair
                long _task_time = _task_queue.front().first;
                _task_queue.pop();
//                _isfree = true;
                _thread_locker.unlock();
//                cout<<"taskID.second: "<<_taskid.second<<endl;
                if (_taskid.second == QUERY) {
//                    cout<<"start querying..."<<endl;
                    long start = clock();

                    vector<KNode> kNNs;
                    if(is_thresholded)
                        kNNs= handleQueryThreshold(_threadid, method_name, k, _taskid.first, mems.dist, mems.visited, mems.q,
                                                                   dijkstra_object_map, threshold, query_id);
                    else
                        kNNs = handleQuery(_threadid, method_name, k, _taskid.first, mems.dist, mems.visited, mems.q,
                                                                   dijkstra_object_map, threshold, query_id);
//                    cout<<"query_id: "<<query_id<<endl;
//                    cout<<"knn size: "<<kNNs.size()<<endl;

//                    vector<KNode> kNNs = DijkstraKNNQuery(k, _taskid.first, mems.dist, mems.visited, mems.q,
//                                                                   dijkstra_object_map);
//                    cout<<"q time: "<<clock()-start<<endl;
//                    cout<<"knn size: "<<kNNs.size()<<endl;
//                    if (kNNs.size() == k && part_maps[_taskid.first]==_threadid) {
                    if (kNNs.size() == k) {
                        global_locker.lock();

                        if(kNNs[k - 1].dis < threshold[query_id])
                            threshold[query_id] = kNNs[k - 1].dis;
//                        cout<<"threshold: "<<threshold[query_id]<<endl;
                        global_locker.unlock();
                    }
                    pair<int, vector<KNode> > partial_res = make_pair(query_id, kNNs);
                    query_id = (query_id + 1) % QUERY_ID_FOLD;
                    _aggregate_thread->add_task(_task_time, partial_res);
//                    cout<<"aggregate task added!"<<endl;

                    // put partial querying result to somewhere
                }
                if (_taskid.second == INSERT) {
                    handleInsert(_threadid, method_name, k, _taskid.first, dijkstra_object_map, mems);

                }
                if (_taskid.second == DELETE) {
                    handleDelete(_threadid, method_name, k, _taskid.first, dijkstra_object_map, mems);

                }

            }
//            else if (_stop) {
            _thread_locker.lock();
            bool cond_stop = _task_queue.empty() && _stop;
            _thread_locker.unlock();
//            if(!cond_stop && needjoin){
//                cout<<"size: "<<_task_queue.size()<<endl;
//            }
            if(cond_stop){
//                cout<<"hi there!!"<<endl;
//                cout<<"final size: "<<_task_queue.size()<<endl;
//                cout << "stopped running" << endl;
                global_locker_2.lock();
                partition_stop_run_threads++;
                global_locker_2.unlock();
//                std::this_thread::sleep_for(std::chrono::microseconds(10));

                break;
            }

        }
//        cout << endl << endl << "out 2" << endl;
    }
};

class PartitionThreadPool {
private:

    std::vector<PartitionThread *> _pool;
    // std::recursive_mutex _locker;
    std::mutex _locker;
    std::thread _main_thread;
    PartitionAggregateThread *_aggregate_thread;
    int num_threads_update;
    int _needjoin;
    double alpha;
    int test_n;
    double fail_p;
    int k;
    int simulation_time;
    vector<int> part_maps;
    double query_rate;
    double insert_rate;
    double delete_rate;
    int num_objects;
    // Algorithm data structure
//    vector<int *> dijkstra_object_map_vec;
public:
    PartitionThreadPool(vector<int> &part_maps_val, int num_threads_update_val, double alpha_val, int k_val,
                        double fail_p_val, int test_n_val, double query_rate_val, double insert_rate_val, double delete_rate_val,
    int simulation_time_val, int num_objects_val) {
//        cout << "constructing PartitionThreadPool..." << endl;
        query_rate = query_rate_val;
        insert_rate = insert_rate_val;
        delete_rate = delete_rate_val;
        part_maps = part_maps_val;
        test_n = test_n_val;
        num_threads_update = num_threads_update_val;
        alpha = alpha_val;
        fail_p = fail_p_val;
        k = k_val;
        simulation_time = simulation_time_val;
        num_objects=num_objects_val;
        threshold.clear();
        for (int i = 0; i <= QUERY_ID_FOLD; i++) threshold.push_back(INT_MAX);
        partition_stop_run_threads=0;

        init();
        int i = num_threads_update;


//        int k_star = compute_k_star(k, num_threads_update, alpha, fail_p);
//        cout<<"k_star: "<<k_star<<endl;
        _aggregate_thread = new PartitionAggregateThread(k, num_threads_update);
//        cout<<"new aggregation thread"<<endl;
        while (i--) {
            PartitionThread *t = new PartitionThread(i, part_maps, k, _aggregate_thread);
            _pool.push_back(t);
        }
        _needjoin = 0;
    }

    //释放线程池
    ~PartitionThreadPool() {
        cout << "hi, PartitionThreadPool()" << endl;
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
        _main_thread.join();
        for (int i = 0; i < _pool.size(); i++) {
            _pool[i]->join();
        }
        _aggregate_thread->join();
    }

    int isNeedJoin() {
        return _needjoin;
    }

    void releasePool() {
        for (int i = 0; i < _pool.size(); ++i) {
            delete _pool[i];
        }
    }

    void start() {

        _main_thread = std::thread(&PartitionThreadPool::run, this);
    }

    void run() {
        mem_struct mems;
        int* car_nodes;
        if(VERIFY){
            allocate_mem(mems, test_n+1);
            car_nodes = new int[test_n+1];
            memset(car_nodes, 0, sizeof(int)*(test_n+1));
        }

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

        for (int i = 0; i < test_n; i++) {

            non_object_list.push_back(i);

        }

//        int simulation_time = 10;
        vector<std::pair<double, int> > full_list;
        for (int i = 0; i < num_objects; i++) {
            full_list.push_back(make_pair(0.0, INSERT));
        }
        vector<std::pair<double, int> > append_list = make_online_query_update_list(query_rate, insert_rate,
                                                                                    delete_rate,
                                                                                    simulation_time);
        for (pair<double, int> &item : append_list)
            full_list.push_back(item);
//        cout << "full_list made..." << endl;
//        cout << "full list size: " << full_list.size() << endl;
        double last_time = 0.0;

        gettimeofday(&global_start, NULL);

        for (int i = 0; i < full_list.size(); i++) {
//            std::this_thread::sleep_for(std::chrono::microseconds(1));
            pair<double, int> &event = full_list[i];
            long issue_time = floor(event.first * MICROSEC_PER_SEC);
            if(!VERIFY) {
                long current_time;
                do {
                    gettimeofday(&end, NULL);
                    current_time = (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
                    if (issue_time <= current_time) break;
//                else std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                } while (true);

//                gettimeofday(&end, NULL);
//                long current_time = (end.tv_sec - global_start.tv_sec) * MICROSEC_PER_SEC + end.tv_usec - global_start.tv_usec;
//                if(current_time < issue_time)
//                    std::this_thread::sleep_for(std::chrono::microseconds(issue_time-current_time));
            }

            // if insert

            if (event.second == INSERT && object_list.size() < test_n) {
//                cout<<"INSERT"<<endl;

                int index = rand() % non_object_list.size();
                int non_object_node = non_object_list[index];
                int last_index = non_object_list.size() - 1;

                int tmp = non_object_list[index];
                non_object_list[index] = non_object_list[last_index];
                non_object_list[last_index] = tmp;
                non_object_list.pop_back();

                object_list.push_back(non_object_node);
                int j_mod = part_maps[non_object_node];
//                cout<<"j_mod: "<<j_mod<<endl;
                // assign insert tasks

                pair<int, int> node_type_pair = std::make_pair(non_object_node, INSERT);
                pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
                _pool[j_mod]->add_task(task);
//                cout<<"insert "<<non_object_node<<endl;
                if(VERIFY){
//                    car_nodes[non_object_node]=1;
                    DijkstraKNNInsert(non_object_node, car_nodes);
                }


            }
            if (event.second == DELETE && object_list.size() > 0) {
//                cout<<"DELETE"<<endl;

                int index = rand() % object_list.size();
                int object_node = object_list[index];
                int last_index = object_list.size() - 1;
                int tmp = object_list[index];
                object_list[index] = object_list[last_index];
                object_list[last_index] = tmp;
                object_list.pop_back();

                non_object_list.push_back(object_node);


                int j_mod = part_maps[object_node];
                pair<int, int> node_type_pair = std::make_pair(object_node, DELETE);
                pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
                _pool[j_mod]->add_task(task);
                if(VERIFY){
//                    car_nodes[object_node]=0;
                    DijkstraKNNDelete(object_node, car_nodes);
                }

            }
//            cout<<"here 2"<<endl;
            if (event.second == QUERY) {
                number_of_queries++;
//                cout<<"QUERY"<<endl;

//                int query_node;
                do {
//                    query_node = rand()%(test_n-1)+1;
                    query_node+=query_gap;
                    query_node = query_node % (test_n-1)+1;

                }while( vStart[query_node] == 0);

                if(VERIFY){
                    vector<KNode> result = DijkstraKNNQuery(k, query_node, mems.dist,
                                                            mems.visited, mems.q, car_nodes);
                    verify_results.push_back(result);

                }
                // put to query tasks

                for (int j = 0; j < num_threads_update; j++) {
                    pair<int, int> node_type_pair = std::make_pair(query_node, QUERY);
                    pair<long, pair<int, int> > task = std::make_pair(issue_time, node_type_pair);
                    _pool[j]->add_task(task);
                }

            }



        }
        if(VERIFY){
            delete[] car_nodes;
            delete_mems(mems);
        }

        for (int i = 0; i < num_threads_update; i++) {
            _pool[i]->set_stop();
        }
        _aggregate_thread->set_stop();


        _needjoin=1;



        cout<<"all set stopped!"<<endl;
    }

};

#endif //TOAIN_PARTITIONQUERYINGTHREAD_H
