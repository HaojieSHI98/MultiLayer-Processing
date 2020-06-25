//
// Created by lostrong on 18/3/10.
//

#ifndef TOAIN_SIMULATION_H
#define TOAIN_SIMULATION_H

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

#endif //TOAIN_SIMULATION_H
