//
// Created by lostrong on 18/5/10.
//

#ifndef TOAIN_MULTICORETRAIN_H
#define TOAIN_MULTICORETRAIN_H


#include<iostream>
#include<vector>
#include <chrono>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include<queue>
#include "../MultiCore/RandomQueryingThread.h"
//#include "MultiCore/PartitionQueryingThread.h"
#include "../data_format/MetisFormat.h"
#include "../scob/ScobConfiguration.h"


// include order is important
#include "../para/GlobalVariables.h"
#include "../RoadKNNUtilities.h"
#include "../dijk/NaiveKNN.h"
#include "../scob/ScobIndex.h"
#include "../util/DataFunctions.h"
#include "../testing/TestToainFreely.h"
#include "../training/ToainTraining.h"
#include "../config/Config.h"
#include "../training/ToainSimulation.h"

int get_num_threads_query(int num_threads_update, int is_single_aggregate, int num_total_threads){
    int num_threads_query = 0;
    if(!is_single_aggregate) {
        num_threads_query = (num_total_threads - 1) / (num_threads_update + 1);// copies
        //   if((multiTestPara.num_total_threads-1)%(num_threads_update+1)) return;
    }
    else {
        num_threads_query = (num_total_threads - 2) / (num_threads_update);// copies
        //   if((multiTestPara.num_total_threads-2)%(num_threads_update)) return;
    }
    cout<<"num_total_threads: "<<num_total_threads<<endl;
    cout<<"num_threads_update: "<<num_threads_update<<endl;
    if(num_threads_query==0) num_threads_query++;
    return num_threads_query;

}

MultiTestParameter anaylizeParameters(){
    MultiTestParameter multiTestPara;

    multiTestPara.num_total_threads = atoi(paras["-multicore"].c_str());

    if(paras.find("-configtime")==paras.end()){
        cout<<"miss configtime parameters!"<<endl;
        stop_here();
    }
    multiTestPara.config_simulation_time = atoi(paras["-configtime"].c_str());// simulation time (seconds)

    if(paras.find("-configstr")==paras.end()){
        cout<<"miss configstr parameters!"<<endl;
        stop_here();
    }
    multiTestPara.configstr = paras["-configstr"].c_str();// simulation time (seconds)
    cout<<"Config String!!!"<<endl;
    cout<<multiTestPara.configstr<<endl<<endl;

    if(paras.find("-testtime")==paras.end()){
        cout<<"miss testtime parameters!"<<endl;
        stop_here();
    }
    multiTestPara.test_simulation_time = atoi(paras["-testtime"].c_str());// simulation time (seconds)

    if(paras.find("-testtime2")==paras.end()){
        multiTestPara.test_simulation_time2 = multiTestPara.test_simulation_time;
    }
    multiTestPara.test_simulation_time2 = atoi(paras["-testtime2"].c_str());// simulation time (seconds)

    if(paras.find("-querythread")==paras.end()){
        multiTestPara.query_thread = multiTestPara.num_total_threads-1;
    }
    else{
        multiTestPara.query_thread = atoi(paras["-querythread"].c_str());
    }

    if(paras.find("-updatethread")==paras.end()){
        multiTestPara.update_thread = 1;
    }
    else
    {
        multiTestPara.update_thread = atoi(paras["-updatethread"].c_str());
    }

    if(paras.find("-DISPLAY")==paras.end()){
        DISPLAY = 0;
    }
    else
    {
        DISPLAY = atoi(paras["-DISPLAY"].c_str());
        cout<<"DISPLAY:"<<DISPLAY<<endl;
    }

    if(paras.find("-querycost")==paras.end()){
        multiTestPara.query_cost = 0;
    }
    else{
        multiTestPara.query_cost = atoi(paras["-querycost"].c_str());
    }

    if(paras.find("-insertcost")==paras.end()){
        multiTestPara.insert_cost = 0;
    }
    else{
        multiTestPara.insert_cost = atoi(paras["-insertcost"].c_str());
    }

    if(paras.find("-deletecost")==paras.end()){
        multiTestPara.delete_cost = 0;
    }
    else{
        multiTestPara.delete_cost = atoi(paras["-deletecost"].c_str());
    }

    if(paras.find("-network")==paras.end()){
        cout<<"here 1"<<endl;
        network_name="NY";
    }
    else{
        network_name=paras["-network"];
    }
    cout<<"network name: "<<network_name<<" para network: "<<paras["-network"]<<endl;

    if(paras.find("-out")==paras.end()){
        cout<<"outfile is missing"<<endl;
        stop_here();
    }
    multiTestPara.suffix = paras["-out"];

    if(paras.find("-method")==paras.end()){
        cout<<"method is missing"<<endl;
        stop_here();
    }

    multiTestPara.method_name=paras["-method"];

    if(paras.find("-threshold")==paras.end()){
        cout<<"threshold is missing"<<endl;
        stop_here();
    }
    multiTestPara.is_thresholded=atoi(paras["-threshold"].c_str());


    if(paras.find("-parmethod")==paras.end()){
        cout<<"parmethod is missing"<<endl;
        stop_here();
    }
    multiTestPara.parmethod = paras["-parmethod"];

    if(paras.find("-single_aggr")==paras.end()){
        cout<<"signle_aggr is missing"<<endl;
        stop_here();
    }
    multiTestPara.is_single_aggregate=atoi(paras["-single_aggr"].c_str());

    if(paras.find("-part")==paras.end()){
        cout<<"miss partition parameters!"<<endl;
        stop_here();
    }

    if(paras.find("-query")==paras.end()){
        cout<<"miss query rate parameters!"<<endl;
        stop_here();
    }

    multiTestPara.query_rate=stod(paras["-query"].c_str());

    if(paras.find("-insert")==paras.end()){
        cout<<"miss insert rate parameters!"<<endl;
        stop_here();
    }

    multiTestPara.insert_rate=stod(paras["-insert"].c_str());

    if(paras.find("-delete")==paras.end()){
        cout<<"miss delete rate parameters!"<<endl;
        stop_here();
    }

    multiTestPara.delete_rate=stod(paras["-delete"].c_str());

    if(paras.find("-query2")==paras.end()){
        multiTestPara.query_rate2 = multiTestPara.query_rate;
    }else multiTestPara.query_rate2=stod(paras["-query2"].c_str());

    if(paras.find("-insert2")==paras.end()){
        multiTestPara.insert_rate2 = multiTestPara.insert_rate;
    }else multiTestPara.insert_rate2=stod(paras["-insert2"].c_str());

    if(paras.find("-delete2")==paras.end()){
        multiTestPara.delete_rate2 = multiTestPara.delete_rate;
    }else multiTestPara.delete_rate2=stod(paras["-delete2"].c_str());

    if(paras.find("-autoconf")==paras.end()){
        multiTestPara.auto_config=0;
    }
    else
        multiTestPara.auto_config= atoi(paras["-autoconf"].c_str());

    if(paras.find("-pfail")==paras.end()){
        multiTestPara.pfail=0.0;
    } else
        multiTestPara.pfail=stod(paras["-pfail"].c_str());

    if(paras.find("-init")==paras.end()){
        multiTestPara.init_objects=1000;
    }
    else
        multiTestPara.init_objects=atoi(paras["-init"].c_str());

    if(paras.find("-layer")==paras.end()){
        multiTestPara.layer = 1;
    }
    else
        multiTestPara.layer=atoi(paras["-layer"].c_str());

    if(paras.find("-toaintype")==paras.end()){
        multiTestPara.toain_type = TOAIN_UPDATE;
    }
    else{
        if(paras["-toaintype"].compare("u")==0){
            multiTestPara.toain_type = TOAIN_UPDATE;
        } else{
            multiTestPara.toain_type = TOAIN_QUERY;
        }
    }

    multiTestPara.num_threads_update = atoi(paras["-part"].c_str());// partitions
    multiTestPara.num_threads_query=get_num_threads_query(multiTestPara.num_threads_update, multiTestPara.is_single_aggregate,
                                                          multiTestPara.num_total_threads);


    return multiTestPara;

}

int chooseSingleTOAINConfiguration(double fail_p, double alpha, int k,
                                   MultiTestParameter& multiTestPara){
    cout<<"start choosing toain configurations for single core..."<<endl;
    prepare_environment();
    level_num = get_level_num_with_dataspace();
    cout<<"level_num: "<<level_num<<endl;
    vector<ScobConfiguration> configurations;
    if(network_name.compare("NY")==0)
        configurations = ScobConfiguration::generate_all_configurations(level_num);
    else
        configurations = ScobConfiguration::generate_all_update_favor_configurations(level_num);
    double query_rate = multiTestPara.query_rate/multiTestPara.layer;
    double insert_rate = multiTestPara.insert_rate;
    double delete_rate = multiTestPara.delete_rate;


    struct timeval start, end;
    cout<<"query, insert, delete: "<<query_rate<<" "<<insert_rate<<" "<<delete_rate<<endl;
    double total_response_time = 0.0;
    query_gap = (int)(test_n/query_rate);
    double smallest_avg_response_time=INT_MAX*1.0;
    int best_conf=-1;

    for(int i = configurations.size()/2;i < configurations.size();i++){
        ScobConfiguration scobConfiguration = configurations[i];
        tradeoff_level = scobConfiguration.tradeoff_level;
        cut_level = scobConfiguration.cut_level;
        mode = scobConfiguration.mode;
        prepare_reverse_shortcut(test_n, tradeoff_level, cut_level, mode);
        cout << "num of threads: " << multiTestPara.num_threads_update << endl;
        multiTestPara.num_threads_query=1;
        multiTestPara.num_threads_update=1;
        multiTestPara.num_total_threads=2;
        cout << "copies: " << multiTestPara.num_threads_query << endl;

        cout<<"testing configuration "<<i<<endl;
        if (multiTestPara.parmethod.compare("rand") == 0) {

            RandomThreadPool *tp = new RandomThreadPool(0, 0, test_n, multiTestPara.num_threads_query, multiTestPara.num_threads_update,
                                                        alpha, k, fail_p, test_n,
                                                        query_rate, insert_rate,
                                                        delete_rate, multiTestPara.config_simulation_time/5,
                                                        multiTestPara.query_cost,multiTestPara.insert_cost,multiTestPara.delete_cost);
            tp->start(); //run the thread

            while (true) {
                std::this_thread::sleep_for(std::chrono::microseconds(1));
                if (tp->isNeedJoin()) {
                    if(can_estimate)
                        gettimeofday(&start, NULL); //start time clock
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&start, NULL);
                        estimate_mutex.unlock();
                    }
                    tp->join();
                    break;
                }
            }
        } else {
        }

        if(can_estimate)
            gettimeofday(&end, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&end, NULL);
            estimate_mutex.unlock();
        }
        cout << end.tv_sec << " " << start.tv_sec << endl;
        cout << "finish in : " << end.tv_sec - start.tv_sec << " secs" << endl;

        if(!multiTestPara.is_single_aggregate) {
            for (int i = 0; i < multiTestPara.num_threads_query; i++) {
                total_response_time += globalThreadVar[0][i]->total_query_time;
                number_of_queries += globalThreadVar[0][i]->number_of_queries;
            }
        }
        else{
            total_response_time += globalThreadVar[0][0]->total_query_time;
            number_of_queries += globalThreadVar[0][0]->number_of_queries;
        }

        cout << "expected response time: " << total_response_time / number_of_queries << " seconds" << endl;
        cout << "total_response_time: " << total_response_time << endl;
        cout << "update response time: "<<total_update_response_time/ number_of_updates<<endl;
        cout << "number_of_queries: " << number_of_queries << endl;
        cout << "query process time: "<<total_query_process_time / number_of_query_processings <<" ";

        cout << "expected_update_response_time: " << total_update_response_time / number_of_updates<< endl;
        cout << "expected_update_process_time: " << total_update_process_time / number_of_updates<< endl;

        if(overload_flag==0 && smallest_avg_response_time > total_response_time / number_of_queries){
            smallest_avg_response_time =  total_response_time / number_of_queries;
            best_conf=i;
        }

        total_response_time = 0.0;
        number_of_queries = 1;
        total_update_response_time = 0.0;
        number_of_updates = 1;
        total_update_process_time=0.0;
        total_query_process_time=0.0;
        number_of_query_processings=1;
    }
    return best_conf;
}
void ChooseTOAINPRConfiguration(double fail_p, double alpha, int k,
                                              int configurationId, MultiTestParameter& multiTestPara){
    cout<<"start choosing configurations..."<<endl;
    is_simulation = 1;
    prepare_environment_with_scob_conf(configurationId);
//    MultiTestParameter multiTestPara = anaylizeParameters();
    struct timeval start, end;
    double query_rate = multiTestPara.query_rate/multiTestPara.layer;
    double insert_rate = multiTestPara.insert_rate;
    double delete_rate = multiTestPara.delete_rate;
    cout<<"query, insert, delete: "<<query_rate<<" "<<insert_rate<<" "<<delete_rate<<endl;
    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;
    double total_response_time = 0.0;
    query_gap = (int)(test_n/query_rate);
    double smallest_avg_response_time=INT_MAX*1.0;
    int best_partition_conf=-1;
    double threshold = ceil(sqrt(multiTestPara.num_total_threads*(insert_rate+delete_rate)/query_rate));
//    if(threshold > multiTestPara.num_total_threads-1) threshold = multiTestPara.num_total_threads-1;
    cout<<"multiTestPara.num_total_threads: "<<multiTestPara.num_total_threads<<endl;
//    cout<<"threshold: "<<threshold<<endl;
//    hier_local_knn_arr_multi = new vector<KNode> *[multiTestPara.num_total_threads];
//    for (int u = 0; u < multiTestPara.num_total_threads; u++) {
//        hier_local_knn_arr_multi[u] = new vector<KNode>[test_n + 1];
//    }
    for(multiTestPara.num_threads_update = 1;multiTestPara.num_threads_update<multiTestPara.num_total_threads &&
            multiTestPara.num_threads_update<=threshold;
        multiTestPara.num_threads_update++) {
        cout << "num of threads: " << multiTestPara.num_threads_update << endl;
        multiTestPara.num_threads_query=get_num_threads_query(multiTestPara.num_threads_update, multiTestPara.is_single_aggregate,
                                                              multiTestPara.num_total_threads);
        cout << "copies: " << multiTestPara.num_threads_query << endl;
//        hier_local_knn_arr_multi = new vector<KNode> *[multiTestPara.num_total_threads];
//        for (int u = 0; u < multiTestPara.num_threads_query; u++) {
////            for (int i = 0; i < multiTestPara.num_threads_update; i++) {
////                hier_local_knn_arr_multi[u * multiTestPara.num_threads_update + i] = new vector<KNode>[test_n + 1];
////            }
//            for (int j = 0; j < multiTestPara.num_threads_update; j++) {
//
//                for (int i = 0; i <= test_n; i++) {
////                    hier_local_knn_arr_multi[u * multiTestPara.num_threads_update + j][i] = hier_local_knn_arr[i];
//                    hier_local_knn_arr_multi[u * multiTestPara.num_threads_update + j][i].clear();
//                }
//            }
//        }

        if (multiTestPara.parmethod.compare("rand") == 0) {

            RandomThreadPool *tp = new RandomThreadPool(0, 0, test_n, multiTestPara.num_threads_query, multiTestPara.num_threads_update,
                                                        alpha, k, fail_p, test_n,
                                                        query_rate, insert_rate,
                                                        delete_rate, multiTestPara.config_simulation_time,
                                                        multiTestPara.query_cost,multiTestPara.insert_cost,multiTestPara.delete_cost);
            tp->start();

            while (true) {
                std::this_thread::sleep_for(std::chrono::microseconds(1));
                if (tp->isNeedJoin()) {
                    if(can_estimate)
                        gettimeofday(&start, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&start, NULL);
                        estimate_mutex.unlock();

                    }
                    tp->join();
                    break;
                }
            }
        } else {
        }

        if(can_estimate)
            gettimeofday(&end, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&end, NULL);
            estimate_mutex.unlock();
        }
        cout << end.tv_sec << " " << start.tv_sec << endl;
        cout << "finish in : " << end.tv_sec - start.tv_sec << " secs" << endl;

        if(!multiTestPara.is_single_aggregate) {
            for (int i = 0; i < multiTestPara.num_threads_query; i++) {
                total_response_time += globalThreadVar[0][i]->total_query_time;
                number_of_queries += globalThreadVar[0][i]->number_of_queries;
            }
        }
        else{
            total_response_time += globalThreadVar[0][0]->total_query_time;
            number_of_queries += globalThreadVar[0][0]->number_of_queries;

        }
        cout << "expected response time: " << total_response_time / number_of_queries << " seconds" << endl;
        cout << "total_response_time: " << total_response_time << endl;
        cout << "update response time: "<<total_update_response_time/ number_of_updates<<endl;
        cout << "number_of_queries: " << number_of_queries << endl;
        cout << "query process time: "<<total_query_process_time / number_of_query_processings <<" ";

        cout << "expected_update_response_time: " << total_update_response_time / number_of_updates<< endl;
        cout << "expected_update_process_time: " << total_update_process_time / number_of_updates<< endl;

        if(overload_flag==0 && smallest_avg_response_time > total_response_time / number_of_queries){
            smallest_avg_response_time =  total_response_time / number_of_queries;
            best_partition_conf=multiTestPara.num_threads_update;
        }
        total_response_time = 0.0;
        number_of_queries = 1;
        total_update_response_time = 0.0;
        number_of_updates = 1;
        total_update_process_time=0.0;
        total_query_process_time=0.0;
        number_of_query_processings=1;
    }
    std::ofstream config_outfile;
    config_outfile.open(input_parameters.output_data_dir + "config_outfile"+(multiTestPara.suffix), std::ios_base::app);
    config_outfile<<"cores: "<< multiTestPara.num_total_threads<<" fail_p: "<<fail_p<<" query: "<<query_rate<<" insert: "<<insert_rate<<" delete: "<<delete_rate
                  <<" best_conf: "<<best_partition_conf<<endl;
    config_outfile.close();


    multiTestPara.num_threads_update=best_partition_conf;
    multiTestPara.num_threads_query=get_num_threads_query(multiTestPara.num_threads_update, multiTestPara.is_single_aggregate,
                                                          multiTestPara.num_total_threads);
//    return multiTestPara;
}

void ChooseVtreePRConfiguration(double fail_p, double alpha, int k, MultiTestParameter& multiTestPara){
    cout<<"start choosing configurations..."<<endl;
    is_simulation = 1;
    struct timeval start, end;

    double query_rate = multiTestPara.query_rate/multiTestPara.layer;
    double insert_rate = multiTestPara.insert_rate;
    double delete_rate = multiTestPara.delete_rate;
    cout<<"query, insert, delete: "<<query_rate<<" "<<insert_rate<<" "<<delete_rate<<endl;
    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;
    double total_response_time = 0.0;
    query_gap = (int)(test_n/query_rate);
    double smallest_avg_response_time=INT_MAX*1.0;
    int best_partition_conf=-1;
    double threshold = ceil(sqrt(multiTestPara.num_total_threads*(insert_rate+delete_rate)/query_rate));
    cout<<"multiTestPara.num_total_threads: "<<multiTestPara.num_total_threads<<endl;
    for(multiTestPara.num_threads_update = 1;multiTestPara.num_threads_update<multiTestPara.num_total_threads &&
                                             multiTestPara.num_threads_update<=threshold;
        multiTestPara.num_threads_update++) {
        memset(car_id_insert_cnt, 0, sizeof(int)*MAX_CORES);

        cout << "num of threads: " << multiTestPara.num_threads_update << endl;
        multiTestPara.num_threads_query=get_num_threads_query(multiTestPara.num_threads_update, multiTestPara.is_single_aggregate,
                                                              multiTestPara.num_total_threads);
        cout << "copies: " << multiTestPara.num_threads_query << endl;

//        string loadfile = input_parameters.input_data_dir + network_name + ".vtree";
        string loadfile = input_parameters.input_data_dir + network_name + ".vtree";
        if(network_name.compare("BJ-old")==0){
            loadfile = input_parameters.input_data_dir  + "BJ.vtree";

        }
        vtrees.clear();
        for (int i = 0; i < multiTestPara.num_threads_query * multiTestPara.num_threads_update; i++) {
            cout << "loading " << i << endl;
//            G_Tree test;
            load_binary(loadfile, i);
//                cout<<"size "<<i<<": "<<sizeof(tree[i])<<endl;

//            vtrees.push_back(&test);
            vtrees.push_back(&(tree[i]));
        }

        if (multiTestPara.parmethod.compare("rand") == 0) {
            RandomThreadPool *tp = new RandomThreadPool(0, 0, test_n, multiTestPara.num_threads_query, multiTestPara.num_threads_update,
                                                        alpha, k, fail_p, test_n,
                                                        query_rate, insert_rate,
                                                        delete_rate, multiTestPara.config_simulation_time,
                                                        multiTestPara.query_cost,multiTestPara.insert_cost,multiTestPara.delete_cost);
            tp->start();

            while (true) {
                std::this_thread::sleep_for(std::chrono::microseconds(1));
                if (tp->isNeedJoin()) {
                    if(can_estimate)
                        gettimeofday(&start, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&start, NULL);
                        estimate_mutex.unlock();

                    }
                    tp->join();
                    break;
                }
            }
        } else {
        }


        if(can_estimate)
            gettimeofday(&end, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&end, NULL);
            estimate_mutex.unlock();
        }
        cout << end.tv_sec << " " << start.tv_sec << endl;
        cout << "finish in : " << end.tv_sec - start.tv_sec << " secs" << endl;

        for (int i = 0; i < multiTestPara.num_threads_query; i++) {
            total_response_time += globalThreadVar[0][i]->total_query_time;
            number_of_queries += globalThreadVar[0][i]->number_of_queries;
        }
        cout << "expected response time: " << total_response_time / number_of_queries << " seconds" << endl;
        cout << "total_response_time: " << total_response_time << endl;
        cout << "update response time: "<<total_update_response_time/ number_of_updates<<endl;
        cout << "number_of_queries: " << number_of_queries << endl;
        cout << "query process time: "<<total_query_process_time / number_of_query_processings <<" ";

        cout << "expected_update_response_time: " << total_update_response_time / number_of_updates<< endl;
        cout << "expected_update_process_time: " << total_update_process_time / number_of_updates<< endl;

        if(overload_flag==0 && smallest_avg_response_time > total_response_time / number_of_queries){
            smallest_avg_response_time =  total_response_time / number_of_queries;
            best_partition_conf=multiTestPara.num_threads_update;
        }
        total_response_time = 0.0;
        number_of_queries = 1;
        total_update_response_time = 0.0;
        number_of_updates = 1;
        total_update_process_time=0.0;
        total_query_process_time=0.0;
        number_of_query_processings=1;

//        for(int i=0;i<multiTestPara.num_threads_query*multiTestPara.num_threads_update;i++){
//            for(int j =0;j<test_n;j++){
//                if(vtree_objects[i][j]){
//                    (vtrees[i])->del_car(j, (vtrees[i])->car_in_node[j][0]);
//                    vtree_objects[i][j]=0;
//                }
//            }
//
//        }
    }
    std::ofstream config_outfile;
    config_outfile.open(input_parameters.output_data_dir + "config_outfile"+(multiTestPara.suffix), std::ios_base::app);
    config_outfile<<"cores: "<< multiTestPara.num_total_threads<<" fail_p: "<<fail_p<<" query: "<<query_rate<<" insert: "<<insert_rate<<" delete: "<<delete_rate
                  <<" best_conf: "<<best_partition_conf<<endl;
    config_outfile.close();


    multiTestPara.num_threads_update=best_partition_conf;
    multiTestPara.num_threads_query=get_num_threads_query(multiTestPara.num_threads_update, multiTestPara.is_single_aggregate,
                                                          multiTestPara.num_total_threads);

}


void ChooseDIJKPRConfiguration(double fail_p, double alpha, int k, MultiTestParameter& multiTestPara){
    cout<<"start choosing configurations..."<<endl;
    srand(time(NULL));
    is_simulation = 1;
    read_road_network();
    struct timeval start, end;
    double query_rate = multiTestPara.query_rate/multiTestPara.layer;
    double insert_rate = multiTestPara.insert_rate;
    double delete_rate = multiTestPara.delete_rate;
    cout<<"query, insert, delete: "<<query_rate<<" "<<insert_rate<<" "<<delete_rate<<endl;
    cout << endl << "tradeoff: " << tradeoff_level << " cut_level: " << cut_level << " mode: " << mode << endl;
    double total_response_time = 0.0;
    query_gap = (int)(test_n/query_rate);
    double smallest_avg_response_time=INT_MAX*1.0;
    int best_partition_conf=-1;
    double threshold = ceil(sqrt(multiTestPara.num_total_threads*(insert_rate+delete_rate)/query_rate));
    cout<<"multiTestPara.num_total_threads: "<<multiTestPara.num_total_threads<<endl;
    for(multiTestPara.num_threads_update = 1;multiTestPara.num_threads_update<multiTestPara.num_total_threads &&
                                             multiTestPara.num_threads_update<=threshold;
        multiTestPara.num_threads_update++) {
        cout << "num of threads: " << multiTestPara.num_threads_update << endl;
        multiTestPara.num_threads_query=get_num_threads_query(multiTestPara.num_threads_update, multiTestPara.is_single_aggregate,
                                                              multiTestPara.num_total_threads);
        cout << "copies: " << multiTestPara.num_threads_query << endl;
        if (multiTestPara.parmethod.compare("rand") == 0) {

            RandomThreadPool *tp = new RandomThreadPool(0, 0, test_n, multiTestPara.num_threads_query, multiTestPara.num_threads_update,
                                                        alpha, k, fail_p, test_n,
                                                        query_rate, insert_rate,
                                                        delete_rate, multiTestPara.config_simulation_time,
                                                        multiTestPara.query_cost,multiTestPara.insert_cost,multiTestPara.delete_cost);
            tp->start();

            while (true) {
                std::this_thread::sleep_for(std::chrono::microseconds(1));
                if (tp->isNeedJoin()) {
                    if(can_estimate)
                        gettimeofday(&start, NULL);
                    else{
                        estimate_mutex.lock();
                        gettimeofday(&start, NULL);
                        estimate_mutex.unlock();

                    }
//                    long remain_time = tp->get_max_remain_time();
                    tp->join();
                    break;
                }
            }
        } else {
        }


        if(can_estimate)
            gettimeofday(&end, NULL);
        else{
            estimate_mutex.lock();
            gettimeofday(&end, NULL);
            estimate_mutex.unlock();
        }
        cout << end.tv_sec << " " << start.tv_sec << endl;
        cout << "finish in : " << end.tv_sec - start.tv_sec << " secs" << endl;

        for (int i = 0; i < multiTestPara.num_threads_query; i++) {
            total_response_time += globalThreadVar[0][i]->total_query_time;
            number_of_queries += globalThreadVar[0][i]->number_of_queries;
        }
        cout << "expected response time: " << total_response_time / number_of_queries << " seconds" << endl;
        cout << "total_response_time: " << total_response_time << endl;
        cout << "number_of_queries: " << number_of_queries << endl;
        cout << "query process time: "<<total_query_process_time / number_of_query_processings <<" ";

        cout << "expected_update_response_time: " << total_update_response_time / number_of_updates<< endl;
        cout << "expected_update_process_time: " << total_update_process_time / number_of_updates<< endl;

        if(overload_flag==0 && smallest_avg_response_time > total_response_time / number_of_queries){
            smallest_avg_response_time =  total_response_time / number_of_queries;
            best_partition_conf=multiTestPara.num_threads_update;
        }
        total_response_time = 0.0;
        number_of_queries = 1;
        total_update_response_time = 0.0;
        number_of_updates = 1;
        total_update_process_time=0.0;
        total_query_process_time=0.0;
        number_of_query_processings=1;
    }
    std::ofstream config_outfile;
    config_outfile.open(input_parameters.output_data_dir + "config_outfile"+(multiTestPara.suffix), std::ios_base::app);
    config_outfile<<"cores: "<< multiTestPara.num_total_threads<<" fail_p: "<<fail_p<<" query: "<<query_rate<<" insert: "<<insert_rate<<" delete: "<<delete_rate
                  <<" best_conf: "<<best_partition_conf<<endl;
    config_outfile.close();


    multiTestPara.num_threads_update=best_partition_conf;
    multiTestPara.num_threads_query=get_num_threads_query(multiTestPara.num_threads_update, multiTestPara.is_single_aggregate,
                                                          multiTestPara.num_total_threads);
//    return multiTestPara;
}
#endif //TOAIN_MULTICORETRAIN_H
