//
// Created by lostrong on 17/7/10.
//

#ifndef SOB_TESTPARAMETERSETTING_H
#define SOB_TESTPARAMETERSETTING_H
#include "../para/GlobalVariables.h"
extern TestParameter test_parameters;
extern InputParameter input_parameters;

void set_test_para() {
    // simulations of command lines

    test_parameters.construct_index = 0;
    test_parameters.verify = 0;
    test_parameters.arrival_model = 0;	//0 FIFO;  1 QF
    test_parameters.data=0; // 2 NW_POI; 1 Beijing_cars
    test_parameters.train_toain_query_first=0;
    test_parameters.train_toain_fifo=0;
    test_parameters.simulate_toain_query_first=0;
    test_parameters.simulate_toain_fifo=0;
    test_parameters.test_toain_query_first=0;
    test_parameters.test_toain_fifo=0;
    test_parameters.test_dijk_fifo=0;
    test_parameters.test_dijk_query_first=0;
    test_parameters.demo=0;
    test_parameters.multi_core=1;
    test_parameters.format_metis=0;

    int count_useful =0;
    if(test_parameters.construct_index) count_useful++;
    if(test_parameters.train_toain_fifo) count_useful++;
    if(test_parameters.train_toain_fifo) count_useful++;
    if(test_parameters.simulate_toain_query_first) count_useful++;
    if(test_parameters.simulate_toain_fifo) count_useful++;
    if(test_parameters.test_dijk_query_first) count_useful++;
    if(test_parameters.test_dijk_fifo) count_useful++;
    if(test_parameters.test_toain_query_first) count_useful++;
    if(test_parameters.test_toain_fifo) count_useful++;
    if(test_parameters.multi_core) count_useful++;
    if(count_useful==0) cout<<"no action is performed!"<<endl;
    if(count_useful>1) cout<<"Too many actions are performed. This may lead to incorrectness"<<endl;


//    input_parameters.set_graph_input("USA-road-d."+network_name+".co", "USA-road-t."+network_name+".gr");
    if(network_name.compare("BJ-old")==0)
        input_parameters.set_graph_input("bj1_node.csv", "bj1_edge.csv");
    else
        input_parameters.set_graph_input("USA-road-d."+network_name+".co", "USA-road-t."+network_name+".gr");


    input_parameters.set_in_out_folder("/home/stone/Downloads/code-MPR/code-MPR/TOAIN/RoadKNN_input/",
                                       "/home/stone/Downloads/code-MPR/code-MPR/TOAIN/RoadKNN_output");
//    input_parameters.set_in_out_folder("/data/sqluo/RoadKNN_input/",
//                                       "/data/sqluo/RoadKNN_output/");
    data_space=data_map[input_parameters.graph_node_file];
    network_name=network_map[input_parameters.graph_node_file];

}
#endif //SOB_TESTPARAMETERSETTING_H
