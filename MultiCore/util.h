//
// Created by lostrong on 18/3/23.
//

#ifndef TOAIN_UTIL_H
#define TOAIN_UTIL_H
#include "../graph/CellNode.h"
#include <iostream>
#include <fstream>

void stop() {
    cout << "input a char to go on ..." << endl;
    char c;
    cin >> c;
}

int verify_two_results_simple(vector<KNode> &r1, vector<KNode> &r2, int k) {

    if (r1.size() != r2.size()) {
        cout << "r1 size != r2 size..." << endl;
        cout<<r1.size()<<" "<<r2.size()<<endl;
        stop();
    }
    for (int i = 0; i < k && i < r1.size(); ++i) {
        if (r1[i].dis != r2[i].dis) {
            int id1, id2;
            id1 = r1[i].id;
            id2 = r2[i].id;

            for (int j = 0; j < r1.size(); ++j) {
                    cout << r1[j].id << " " << r1[j].dis
                         << endl;
                    cout << r2[j].id << " " << r2[j].dis << endl << endl;

            }
            stop();
        }

    }
    return 1;
}

void gen_random_numbers(){
    cout<<"generating random numbers..."<<endl;
        fstream _file;
        _file.open(input_parameters.output_data_dir + "random_numbers", ios::in);
        int file_exist=1;
        if(!_file) {
            file_exist = 0;
        }
        if(!file_exist){
            cout<<"file not exists"<<endl;
            _file.close();
            _file.open(input_parameters.output_data_dir + "random_numbers", ios::out);
            for(int i = 0;i<100000000;i++){
                int j=rand();
                global_random_numbers[i]=j;
                _file<<j<<endl;

            }
            _file.close();
        }
        else{
            cout<<"file exists"<<endl;

            for(int i = 0;i<100000000;i++){
                _file>>global_random_numbers[i];
            }
            _file.close();

        }

    cout<<"finish generating random numbers..."<<endl;
};
#endif //TOAIN_UTIL_H
