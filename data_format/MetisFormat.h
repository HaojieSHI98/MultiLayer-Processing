//
// Created by lostrong on 18/3/19.
//

#ifndef TOAIN_METISFORMAT_H
#define TOAIN_METISFORMAT_H
#include "../para/GlobalVariables.h"
#include "../io/Input.h"

void transformToMetis(int* vStart, int* eNext, int test_n, int test_m){
    read_road_network();
    string out_path = "BJ_";
    out_path = out_path.append("metis_graph.txt");
    out_path = input_parameters.output_data_dir+out_path;
    FILE* out_file = fopen(out_path.c_str(), "w");
    fprintf(out_file, "%d %d\n", test_n, test_m/2);
    for(int u = 1;u<=test_n;u++){
        vector<int> tmp_list;
        for (int i = vStart[u]; i != 0; i = eNext[i]) {
            int v = eNode[i];
            tmp_list.push_back(v);
        }
        for(int i=0;i<tmp_list.size()-1;i++){
            fprintf(out_file, "%d ", tmp_list[i]);
        }
        if(tmp_list.size()){
            fprintf(out_file, "%d", tmp_list[tmp_list.size()-1]);
        }

        fprintf(out_file, "\n");
    }
    fclose(out_file);


}

#endif //TOAIN_METISFORMAT_H
