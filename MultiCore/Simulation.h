//
// Created by lostrong on 18/3/10.
//

#ifndef TOAIN_SIMULATION_H
#define TOAIN_SIMULATION_H

vector<string> split(const string& str, const string& delim) {
    vector<string> res;
    if("" == str) return res;
    //先将要切割的字符串从string类型转换为char*类型
    char * strs = new char[str.length() + 1] ; //不要忘了
    strcpy(strs, str.c_str());

    char * d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());

    char *p = strtok(strs, d);
    while(p) {
        string s = p; //分割得到的字符串转换为string类型
        res.push_back(s); //存入结果数组
        p = strtok(NULL, d);
    }
    return res;
}
vector<std::pair<double, int> > make_online_query_update_list(double query_rate, double insert_rate, double delete_rate, int simulation_time){
    vector<double> query_sequence = poisson(query_rate, simulation_time);
    vector<double> insert_sequence = poisson(insert_rate, simulation_time);
    vector<double> delete_sequence = poisson(delete_rate, simulation_time);
    int query_size = query_sequence.size();
    int insert_size = insert_sequence.size();
    int delete_size = delete_sequence.size();
    query_sequence.push_back(simulation_time+1);
    insert_sequence.push_back(simulation_time+1);
    delete_sequence.push_back(simulation_time+1);

    int query_index=0, insert_index=0, delete_index=0;

    vector<std::pair<double, int> > full_list;
    while(query_index<query_size|| insert_index<insert_size || delete_index<delete_size){
        if(query_sequence[query_index]<= insert_sequence[insert_index] &&
                query_sequence[query_index] <= delete_sequence[delete_index]){
            full_list.push_back(make_pair(query_sequence[query_index], QUERY));
            query_index++;

        }
        if(insert_sequence[insert_index] <= query_sequence[query_index] &&
                insert_sequence[insert_index] <= delete_sequence[delete_index]){
            full_list.push_back(make_pair(insert_sequence[insert_index], INSERT));
            insert_index++;
        }
        if(delete_sequence[delete_index] <= query_sequence[query_index] &&
           delete_sequence[delete_index] <= insert_sequence[insert_index]){
            full_list.push_back(make_pair(delete_sequence[delete_index], DELETE));
            delete_index++;

        }

    }
    return full_list;


}

vector<std::pair<double, int> > make_online_query_update_list_new(double query_rate1, double insert_rate1, double delete_rate1, int simulation_time1,
                                                                  double query_rate2, double insert_rate2, double delete_rate2, int simulation_time2){
    vector<double> query_sequence = poisson(query_rate1, simulation_time1);
    vector<double> insert_sequence = poisson(insert_rate1, simulation_time1);
    vector<double> delete_sequence = poisson(delete_rate1, simulation_time1);

//    query_sequence.push_back(simulation_time1+1);
//    insert_sequence.push_back(simulation_time1+1);
//    delete_sequence.push_back(simulation_time1+1);
    vector<double> query_sequence2 = poisson(query_rate2, simulation_time2);
    vector<double> insert_sequence2 = poisson(insert_rate2, simulation_time2);
    vector<double> delete_sequence2 = poisson(delete_rate2, simulation_time2);
    for(int j=0;j<query_sequence2.size();j++)
    {
        query_sequence2[j] += simulation_time1;
    }
    for(int j=0;j<insert_sequence2.size();j++)
    {
        insert_sequence2[j] += simulation_time1;
    }
    for(int j=0;j<delete_sequence2.size();j++)
    {
        delete_sequence2[j] += simulation_time1;
    }

    query_sequence.insert(query_sequence.end(),query_sequence2.begin(),query_sequence2.end());
    insert_sequence.insert(insert_sequence.end(),insert_sequence2.begin(),insert_sequence2.end());
    delete_sequence.insert(delete_sequence.end(),delete_sequence2.begin(),delete_sequence2.end());

    int query_size = query_sequence.size();
    int insert_size = insert_sequence.size();
    int delete_size = delete_sequence.size();

    query_sequence.push_back(simulation_time1+simulation_time2+1);
    insert_sequence.push_back(simulation_time1+simulation_time2+1);
    delete_sequence.push_back(simulation_time1+simulation_time2+1);

    int query_index=0, insert_index=0, delete_index=0;
    vector<std::pair<double, int> > full_list;
    while(query_index<query_size|| insert_index<insert_size || delete_index<delete_size){
        if(query_sequence[query_index]<= insert_sequence[insert_index] &&
           query_sequence[query_index] <= delete_sequence[delete_index]){
            full_list.push_back(make_pair(query_sequence[query_index], QUERY));
            query_index++;

        }
        if(insert_sequence[insert_index] <= query_sequence[query_index] &&
           insert_sequence[insert_index] <= delete_sequence[delete_index]){
            full_list.push_back(make_pair(insert_sequence[insert_index], INSERT));
            insert_index++;
        }
        if(delete_sequence[delete_index] <= query_sequence[query_index] &&
           delete_sequence[delete_index] <= insert_sequence[insert_index]){
            full_list.push_back(make_pair(delete_sequence[delete_index], DELETE));
            delete_index++;

        }

    }
    return full_list;


}

vector<std::pair<double, int> > make_online_query_update_list_str(string configstr, int objectnum,int layer){
    vector<string> configstring = split(configstr,"_");
    cout<<"config_String!!"<<endl;
    if(configstring.size()%3!=0) stop_here();
    int test_time_all = 0;
    vector<double> query_sequence,insert_sequence,delete_sequence;
    for(int j_c = 0;j_c<configstring.size();j_c+=3)
    {
        int query_ratio = atoi(configstring[j_c].c_str());
        int update_ratio = atoi(configstring[j_c+1].c_str());
        int test_time = atoi(configstring[j_c+2].c_str());
        int updatenum = objectnum*update_ratio;
        int insertnum = updatenum/2;
        int deletenum = updatenum/2;
        int querynum = objectnum*query_ratio/layer;
        cout<<"update_num: "<<updatenum<<" query_num:"<<querynum<<endl;
        vector<double> query_sequence_sub = poisson(querynum, test_time);
        vector<double> insert_sequence_sub = poisson(insertnum, test_time);
        vector<double> delete_sequence_sub = poisson(deletenum, test_time);
        if(j_c!=0){
            for(int j=0;j<query_sequence_sub.size();j++)
            {
                query_sequence_sub[j] += test_time_all;
            }
            for(int j=0;j<insert_sequence_sub.size();j++)
            {
                insert_sequence_sub[j] += test_time_all;
            }
            for(int j=0;j<delete_sequence_sub.size();j++)
            {
                delete_sequence_sub[j] += test_time_all;
            }
        }
        test_time_all +=test_time;
        query_sequence.insert(query_sequence.end(),query_sequence_sub.begin(),query_sequence_sub.end());
        insert_sequence.insert(insert_sequence.end(),insert_sequence_sub.begin(),insert_sequence_sub.end());
        delete_sequence.insert(delete_sequence.end(),delete_sequence_sub.begin(),delete_sequence_sub.end());
//        cout<<configstring[j_c]<<endl;
        cout<<"query: "<<querynum<<" insert: "<<insertnum<<" delete: "<<deletenum<<" test_time: "<<test_time<<endl;
    }
//    cout<<endl<<endl;
//    stop_here();

    int query_size = query_sequence.size();
    int insert_size = insert_sequence.size();
    int delete_size = delete_sequence.size();

    query_sequence.push_back(test_time_all+1);
    insert_sequence.push_back(test_time_all+1);
    delete_sequence.push_back(test_time_all+1);

    int query_index=0, insert_index=0, delete_index=0;
    vector<std::pair<double, int> > full_list;
    while(query_index<query_size|| insert_index<insert_size || delete_index<delete_size){
        if(query_sequence[query_index]<= insert_sequence[insert_index] &&
           query_sequence[query_index] <= delete_sequence[delete_index]){
            full_list.push_back(make_pair(query_sequence[query_index], QUERY));
            query_index++;

        }
        if(insert_sequence[insert_index] <= query_sequence[query_index] &&
           insert_sequence[insert_index] <= delete_sequence[delete_index]){
            full_list.push_back(make_pair(insert_sequence[insert_index], INSERT));
            insert_index++;
        }
        if(delete_sequence[delete_index] <= query_sequence[query_index] &&
           delete_sequence[delete_index] <= insert_sequence[insert_index]){
            full_list.push_back(make_pair(delete_sequence[delete_index], DELETE));
            delete_index++;

        }

    }
    return full_list;


}

#endif //TOAIN_SIMULATION_H
