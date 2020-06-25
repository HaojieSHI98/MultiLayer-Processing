//
// Created by lostrong on 17/7/10.
//

#ifndef SOB_TESTDEPRECATED_H
#define SOB_TESTDEPRECATED_H
#include <string.h>
#include <iostream>
#include "../para/GlobalVariables.h"
#include "../util/GeneralUtility.h"
#include "../util/TimeStat.h"
#include "../util/DataFunctions.h"
#include "../queue/GeneralQueueSimulator.h"



// deprecated
void testing_poke_one_instance_compare(string method_name, int T_r, double obj_ratio, int k, double update_obj_ratio,
                                       vector<int> &query_costs,
                                       vector<int> &update_costs) {
    simulate_times costs;

    cout << "---------------------------------------------------------" << endl;
    cout << "poke's model:" << endl;


    cout << "T_r is " << T_r << endl;
    cout << "Total objects: " << test_n * obj_ratio << endl;
    cout << "Update arrival rate: " << test_n * obj_ratio / update_obj_ratio << " ";
//    string method_name = "GTree";
    string query_simulate_input_file_name, update_simulate_input_file_name;
    if (test_parameters.data == 0) {
        if (k == 20) {
            query_simulate_input_file_name = make_out_file_name(method_name, -1, obj_ratio, "query_simulate");

        } else {
            query_simulate_input_file_name = make_out_file_name(method_name, k, -1, "query_simulate");

        }
        update_simulate_input_file_name = make_out_file_name(method_name, -1, obj_ratio, "update_simulate");
    } else {
        query_simulate_input_file_name = make_out_file_name(method_name, k, -1, "query_simulate");
        update_simulate_input_file_name = make_out_file_name(method_name, -1, -1, "update_simulate");
    }
    FILE *query_simulate_intput_file = fopen((output_data_dir + query_simulate_input_file_name).c_str(), "r");
    if (query_simulate_intput_file == NULL) {
        cout << query_simulate_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    FILE *update_simulate_input_file = fopen((output_data_dir + update_simulate_input_file_name).c_str(), "r");
    if (update_simulate_input_file == NULL) {
        cout << update_simulate_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    cout << "open query simulating file: " << query_simulate_input_file_name << endl;
    cout << "open update simulating file: " << update_simulate_input_file_name << endl;

    int query_val;
    int update_val;
    while (fscanf(query_simulate_intput_file, "%d", &query_val) == 1)
        query_costs.push_back(query_val);
    fclose(query_simulate_intput_file);
    while (fscanf(update_simulate_input_file, "%d", &update_val) == 1)
        update_costs.push_back(update_val);
    fclose(update_simulate_input_file);

}

// deprecated
void testing_poke_compare(string method_name) {
    test_parameters.arrival_model = 0;
    int vary_k = 0;
    int vary_o = 1;
    int vary_update_ratio = 0;
    int vary_Tr = 0;
    string result_out_file_name;
    if (vary_k == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyk_result");
    if (vary_o == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyo_result");
    if (vary_update_ratio == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyUr_result");
    if (vary_Tr == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyTr_result");
    FILE *result_out_file = fopen((output_data_dir + result_out_file_name).c_str(), "w");
    if (result_out_file == NULL) {
        cout << "cannot open file: " << result_out_file_name << " in test_toain_fifo()" << endl;
        stop_here();
    }
    fprintf(result_out_file, "T_r \t update_arrival \t k \t \t Obj_ratio \t Thr\n");

    vector<double> obj_ratio_vals = {0.0001, 0.001, 0.01, 0.1};
    vector<double> obj_ratio_vals_def = {0.001};
    vector<int> k_vals = {1, 10, 20, 30, 40};
    vector<int> k_vals_def = {20};
    vector<int> Tr_vals = {100, 200, 400, 800, 1600, 3200, 6400};
    vector<int> Tr_vals_def = {800};
    vector<double> Ur_vals = {1.0, 2.0, 4.0, 8.0, 16.0};
    vector<double> Ur_vals_def = {4.0};

    vector<double> selected_objs = obj_ratio_vals_def, selected_Ur_vals = Ur_vals_def;
    vector<int> selected_k_vals = k_vals_def, selected_Tr_vals = Tr_vals_def;
    if (vary_o) selected_objs = obj_ratio_vals;
    if (vary_k) selected_k_vals = k_vals;
    if (vary_update_ratio) selected_Ur_vals = Ur_vals;
    if (vary_Tr) selected_Tr_vals = Tr_vals;

    for (double obj_ratio:selected_objs) {
        for (int k:selected_k_vals) {
            for (int T_r : selected_Tr_vals) {
                for (double update_obj_ratio: selected_Ur_vals) {
                    if (test_parameters.data == 1) obj_ratio = 2770 * 1.0 / test_n;
                    vector<int> query_costs;
                    vector<int> update_costs;
                    testing_poke_one_instance_compare(method_name, T_r, obj_ratio, k, update_obj_ratio, query_costs,
                                                      update_costs);

                    vector<double> simulated_query_costs;
                    vector<double> simulated_update_costs;
                    int th = test_for_fifo(simulated_query_costs, simulated_update_costs,
                                           10, T_r, test_n * obj_ratio / update_obj_ratio, query_costs,
                                           update_costs);
                    fprintf(result_out_file, "%d \t %lf \t %d \t %lf \t %d\n", T_r, update_obj_ratio, k, obj_ratio, th);

                    cout << T_r << " " << update_obj_ratio << " " << k << " " << obj_ratio << " " << th << endl;

                }
            }
        }
    }
    fclose(result_out_file);
}


// deprecated
void testing_uber_compare(string method_name) {
    test_parameters.arrival_model = 1;
    int vary_k = 0;
    int vary_o = 0;
    int vary_T = 1;
    int vary_Tr = 0;
    string result_out_file_name;
    if (vary_k == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyk_result");
    if (vary_o == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyo_result");
    if (vary_T == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyT_result");
    if (vary_Tr == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyTr_result");
    FILE *result_out_file = fopen((output_data_dir + result_out_file_name).c_str(), "w");
    if (result_out_file == NULL) {
        cout << "cannot open file: " << result_out_file_name << " in test_toain_fifo()" << endl;
        stop_here();
    }
    fprintf(result_out_file, "T_r \t update_arrival \t k \t \t Obj_ratio \t Thr\n");

    vector<double> obj_ratio_vals = {0.0001, 0.001, 0.01, 0.1};
    vector<double> obj_ratio_vals_def = {0.001};
    vector<int> k_vals = {1, 10, 20, 30, 40};
    vector<int> k_vals_def = {20};
    vector<int> Tr_vals = {100, 200, 400, 800, 1600, 3200, 6400};
    vector<int> Tr_vals_def = {800};
    vector<double> T_vals = {1.0, 2.0, 4.0, 8.0, 16.0};
    vector<double> T_vals_def = {4.0};

    vector<double> selected_objs = obj_ratio_vals_def, selected_T_vals = T_vals_def;
    vector<int> selected_k_vals = k_vals_def, selected_Tr_vals = Tr_vals_def;
    if (vary_o) selected_objs = obj_ratio_vals;
    if (vary_k) selected_k_vals = k_vals;
    if (vary_T) selected_T_vals = T_vals;
    if (vary_Tr) selected_Tr_vals = Tr_vals;

    for (double obj_ratio:selected_objs) {
        for (int k:selected_k_vals) {
            for (int T_r : selected_Tr_vals) {
                for (double T: selected_T_vals) {
                    if (test_parameters.data == 1) obj_ratio = 2770 * 1.0 / test_n;
                    vector<int> query_costs;
                    vector<int> update_costs;
                    testing_uber_one_instance_compare(method_name, T_r, obj_ratio, k, obj_ratio, query_costs,
                                                      update_costs);

                    vector<double> simulated_query_costs;
                    vector<double> simulated_update_costs;
                    int th = test_for_queryfirst(simulated_query_costs, simulated_update_costs,
                                                 T, T_r, T * MICROSEC_PER_SEC, query_costs, update_costs);

                    fprintf(result_out_file, "%d \t %lf \t %d \t %lf \t %d\n", T_r, obj_ratio, k, obj_ratio, th);

                    cout << T_r << " " << obj_ratio << " " << k << " " << obj_ratio << " " << th << endl;

                }
            }
        }
    }
    fclose(result_out_file);
}


// deprecated
void testing_uber_one_instance_compare(string method_name, int T_r, double obj_ratio, int k, double update_obj_ratio,
                                       vector<int> &query_costs,
                                       vector<int> &update_costs) {
    simulate_times costs;

    cout << "---------------------------------------------------------" << endl;
    cout << "uber's model:" << endl;


    cout << "T_r is " << T_r << endl;
    cout << "Total objects: " << test_n * obj_ratio << endl;
    cout << "Update arrival rate: " << test_n * obj_ratio / update_obj_ratio << " ";
//    string method_name = "GTree";
    string query_simulate_input_file_name, update_simulate_input_file_name;
    if (test_parameters.data == 0) {
        if (k == 20) {
            query_simulate_input_file_name = make_out_file_name(method_name, -1, obj_ratio, "query_simulate");

        } else {
            query_simulate_input_file_name = make_out_file_name(method_name, k, -1, "query_simulate");

        }
        update_simulate_input_file_name = make_out_file_name(method_name, -1, obj_ratio, "update_simulate");
    } else {
        query_simulate_input_file_name = make_out_file_name(method_name, k, -1, "query_simulate");
        update_simulate_input_file_name = make_out_file_name(method_name, -1, -1, "update_simulate");
    }
    FILE *query_simulate_intput_file = fopen((output_data_dir + query_simulate_input_file_name).c_str(), "r");
    if (query_simulate_intput_file == NULL) {
        cout << query_simulate_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    FILE *update_simulate_input_file = fopen((output_data_dir + update_simulate_input_file_name).c_str(), "r");
    if (update_simulate_input_file == NULL) {
        cout << update_simulate_input_file_name << " cannot open!" << endl;
        stop_here();
    }
    cout << "open query simulating file: " << query_simulate_input_file_name << endl;
    cout << "open update simulating file: " << update_simulate_input_file_name << endl;

    int query_val;
    int update_val;
    while (fscanf(query_simulate_intput_file, "%d", &query_val) == 1)
        query_costs.push_back(query_val);
    fclose(query_simulate_intput_file);
    while (fscanf(update_simulate_input_file, "%d", &update_val) == 1)
        update_costs.push_back(update_val);
    fclose(update_simulate_input_file);

}



void testing_dijk_one_instance(string method_name, int T_r, double obj_ratio, int k, double update_obj_ratio,
                               vector<int> &query_costs,
                               vector<int> &update_costs) {
    simulate_times costs;

    cout << "---------------------------------------------------------" << endl;
    cout << "poke's model:" << endl;


    cout << "T_r is " << T_r << endl;
    cout << "Total objects: " << test_n * obj_ratio << endl;
    cout << "Update arrival rate: " << test_n * obj_ratio / update_obj_ratio << " ";
//    string method_name = "GTree";
    string query_simulate_input_file_name;
    if (test_parameters.data == 0) {
        query_simulate_input_file_name = make_out_file_name(method_name, k, obj_ratio, "simulate");


    } else {
        query_simulate_input_file_name = make_out_file_name(method_name, k, -1, "simulate");

    }
    FILE *query_simulate_intput_file = fopen((output_data_dir + query_simulate_input_file_name).c_str(), "r");
    if (query_simulate_intput_file == NULL) {
        cout << query_simulate_input_file_name << " cannot open!" << endl;
        stop_here();
    }

    cout << "open query simulating file: " << query_simulate_input_file_name << endl;

    int query_val;
    while (fscanf(query_simulate_intput_file, "%d", &query_val) == 1)
        query_costs.push_back(query_val);
    fclose(query_simulate_intput_file);


}

// deprecated
void testing_dijk() {
    test_parameters.arrival_model = 0;
    string method_name = "dijk";
    int vary_k = 0;
    int vary_o = 0;
    int vary_update_ratio = 0;
    int vary_Tr = 0;
    string result_out_file_name;
    if (vary_k == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyk_result");
    if (vary_o == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyo_result");
    if (vary_update_ratio == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyUr_result");
    if (vary_Tr == 1)
        result_out_file_name = make_out_file_name(method_name, -1, -1.0, "varyTr_result");
    if (vary_k == 0 && vary_o == 0 && vary_update_ratio == 0 && vary_Tr == 0)
        result_out_file_name = "case_study";
    FILE *result_out_file = fopen((output_data_dir + result_out_file_name).c_str(), "w");
    if (result_out_file == NULL) {
        cout << "cannot open file: " << result_out_file_name << " in test_toain_fifo()" << endl;
        stop_here();
    }
    fprintf(result_out_file, "T_r \t update_arrival \t k \t \t Obj_ratio \t Thr\n");

    vector<double> obj_ratio_vals = {0.001, 0.01, 0.1};
    vector<double> obj_ratio_vals_def = {200000.0 / test_n};
    vector<int> k_vals = {1, 10, 20, 30, 40};
    vector<int> k_vals_def = {20};
    vector<int> Tr_vals = {100, 200, 400, 800, 1600, 3200, 6400};
    vector<int> Tr_vals_def = {800};
    vector<double> Ur_vals = {1.0, 2.0, 4.0, 8.0, 16.0};
    vector<double> Ur_vals_def = {4.0};

    vector<double> selected_objs = obj_ratio_vals_def, selected_Ur_vals = Ur_vals_def;
    vector<int> selected_k_vals = k_vals_def, selected_Tr_vals = Tr_vals_def;
    if (vary_o) selected_objs = obj_ratio_vals;
    if (vary_k) selected_k_vals = k_vals;
    if (vary_update_ratio) selected_Ur_vals = Ur_vals;
    if (vary_Tr) selected_Tr_vals = Tr_vals;

    for (double obj_ratio:selected_objs) {
        for (int k:selected_k_vals) {
            for (int T_r : selected_Tr_vals) {
                for (double update_obj_ratio: selected_Ur_vals) {
                    if (test_parameters.data == 1) obj_ratio = 2770 * 1.0 / test_n;

                    vector<int> query_costs, update_costs;
                    for (int kk = 0; kk < 100; kk++)
                        update_costs.push_back(0);
                    testing_dijk_one_instance(method_name, T_r, obj_ratio, k, update_obj_ratio, query_costs,
                                              update_costs);
                    cout << query_costs.size() << endl;
                    vector<double> simulated_query_costs;
                    vector<double> simulated_update_costs;
                    int th = test_for_fifo(simulated_query_costs, simulated_update_costs,
                                           10, T_r, test_n * obj_ratio / update_obj_ratio, query_costs, update_costs);
                    fprintf(result_out_file, "%d \t %lf \t %d \t %lf \t %d\n", T_r, update_obj_ratio, k, obj_ratio, th);

                    cout << T_r << " " << update_obj_ratio << " " << k << " " << obj_ratio << " " << th << endl;

                }
            }
        }
    }
    fclose(result_out_file);
}

#endif //SOB_TESTDEPRECATED_H
