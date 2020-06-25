/*
 * TestPara.h
 *
 *  Created on: 2017年2月14日
 *      Author: lostrong
 */

#ifndef SRC_TESTPARA_H_
#define SRC_TESTPARA_H_

struct MultiTestParameter{
    int layer;
	int init_objects;
	int num_total_threads;
	int num_threads_update;
	int config_simulation_time;// simulation time (seconds)
    int test_simulation_time;
	string suffix;
	string method_name;
	int is_thresholded;
	string parmethod;
	int	is_single_aggregate;
	int num_threads_query;
	double query_rate;
	double insert_rate;
	double delete_rate;
    int auto_config;
    double pfail;
    int toain_type;
    int query_cost;
    int insert_cost;
    int delete_cost;
    int query_thread;
    int update_thread;
};

struct TestParameter{
	int version;
	int construct_index;
	int verify;
	int arrival_model;
	int data; // 1: shenzhou 0: others
	int train_toain_query_first;
	int train_toain_fifo;
	int simulate_toain_query_first;
	int simulate_toain_fifo;
	int test_toain_query_first;
	int test_toain_fifo;
	int test_dijk_query_first;
	int test_dijk_fifo;
	int demo;
	int multi_core;
	int format_metis;

};

struct InputParameter{
    // we separate the graph files according to the file formats stated in http://www.dis.uniroma1.it/challenge9/download.shtml
    string graph_node_file;
    string graph_edge_file;
    string input_data_dir;
    string output_data_dir;


    void set_graph_input(string graph_node_file, string graph_edge_file){
        this->graph_node_file = graph_node_file;
        this->graph_edge_file = graph_edge_file;

    }

    void set_in_out_folder(string input_data_dir, string out_data_dir){
        this->input_data_dir=input_data_dir;
        this->output_data_dir=out_data_dir;
    }
};



#endif /* SRC_TESTPARA_H_ */
