//
// Created by lostrong on 18/3/23.
//

#ifndef TOAIN_FUNCTIONTEMPLATE_H
#define TOAIN_FUNCTIONTEMPLATE_H


#include "../dijk/DijkstraQueue.h"
#include "../graph/CellNode.h"
#include "../../VTree-master/vtree_knn_demo.h"
#include "DijkstraKNN.h"
#include "../scob/ScobQuery.h"
#include "../scob/ScobIndex.h"
#include <string>


int k_max_for_toain = 50;
int durationTime_for_toain = 5;
//static int car_id_remove_cnt=0;

vector<KNode> HandleQueryThreshold(int copy_id, int thread_id, string method_name, int k, int query_node, long *dist,
                                   int *visited, DijkstraQueue *q,
                                   int *object_map, vector<int> &threshold, int query_id){
    if(method_name.compare("dijk")==0){
        return DijkstraKNNQueryThreshold(k, query_node, dist,
                visited, q, object_map, threshold, query_id);
    }
    if(method_name.compare("vtree")==0){
//        cout<<"handle query threshold"<<endl;
        return vtrees[copy_id*multiTestPara.num_threads_update+thread_id]->KNN_min_dist_car_pair_threshold(query_node, k, threshold, query_id);
//        return vtrees[thread_id]->KNN_min_dist_car_pair(query_node, k);
//        return tree.KNN_min_dist_car_pair(query_node, k);
    }
    if(method_name.compare("toain")==0){
//        cout<<"handle toain query"<<endl;
//        long start = clock();
//        vector<KNode> res = scob_KNN_query(k, hier_local_knn_arr_multi[thread_id], query_node, dist, visited, q);
        vector<KNode> res = scob_KNN_query_threshold(k, hier_local_knn_arr_multi[copy_id*multiTestPara.num_threads_update+thread_id], query_node, dist, visited, q, threshold, query_id);
//        cout<<"time: "<<clock()-start<<endl;
        return res;
    }
}

vector<KNode> HandleQuery(int copy_id, int thread_id, string method_name, int k, int query_node, long *dist,
                          int *visited, DijkstraQueue *q,
                          int *object_map, vector<int> &threshold, int query_id,
                          vector<KNode> *hier_local_knn_arr_single){
    if(method_name.compare("dijk")==0){
        return DijkstraKNNQuery(k, query_node, dist,
                                         visited, q, object_map);
    }
    if(method_name.compare("vtree")==0){
//        cout<<"handle query"<<endl;
//        cout<<"using index copy: "<<copy_id*num_threads_update+thread_id<<endl;
        return vtrees[copy_id*multiTestPara.num_threads_update+thread_id]->KNN_min_dist_car_pair(query_node, k);
//        return vtree.KNN_min_dist_car_pair(query_node, k);
//        return tree.KNN_min_dist_car_pair(query_node, k);
    }
    if(method_name.compare("toain")==0){
//        cout<<"handle toain query"<<endl;
//        cout<<"query node: "<<query_node<<endl;
        return scob_KNN_query(k, hier_local_knn_arr_single, query_node, dist, visited, q);
    }

}

void HandleDelete(int copy_id, int thread_id, string method_name, int k, int object_node,
                  int *object_map, mem_struct &mems, vector<KNode> *hier_local_knn_arr_single) {
    if(method_name.compare("dijk")==0){
        DijkstraKNNDelete(object_node, object_map);
    }
    if(method_name.compare("vtree")==0){
        // always delete the first car in the node
//        cout<<"handle delete"<<endl;
//        cout<<"using index copy: "<<copy_id*num_threads_update+thread_id<<endl;
        (vtrees[copy_id*multiTestPara.num_threads_update+thread_id])->del_car(object_node, (vtrees[copy_id*multiTestPara.num_threads_update+thread_id])->car_in_node[object_node][0]);
//        vtree.del_car(object_node, vtree.car_in_node[object_node][0]);
//        tree.del_car(object_node, (tree).car_in_node[object_node][0]);
        object_map[object_node]=0;
    }
    if(method_name.compare("toain")==0){
//        cout<<"handle toain delete"<<endl;
        remove_kdnn_object_multi(object_node, k_max_for_toain, k, mems, rev_shortcuts_arr, durationTime_for_toain,
                                 hier_local_knn_arr_single);
    }


}

void HandleInsert(int copy_id, int thread_id, string method_name, int k, int object_node,
                  int *object_map, mem_struct &mems, vector<KNode> *hier_local_knn_arr_single){


    if(method_name.compare("dijk")==0){
        DijkstraKNNInsert(object_node, object_map);
    }

    if(method_name.compare("vtree")==0){
//        cout<<"handle vtree insert "<<thread_id<<endl;

//        cout<<"insert using index copy: "<<copy_id*multiTestPara.num_threads_update+thread_id<<endl;
        vtrees[copy_id*multiTestPara.num_threads_update+thread_id]->add_car(object_node, car_id_insert_cnt[copy_id]);
//        vtree.add_car(object_node, car_id_insert_cnt[copy_id]);
//        tree[thread_id].add_car(object_node, car_id_insert_cnt);
//        tree.add_car(object_node, car_id_insert_cnt);
        car_id_insert_cnt[copy_id]=(car_id_insert_cnt[copy_id]+1)%test_n;
//        cout<<"finish vtree insert "<<thread_id<<endl;
        object_map[object_node]=1;
    }
    if(method_name.compare("toain")==0){
//        long start = clock();
//        cout<<"handle toain insert"<<endl;
        if(number_of_updates<=multiTestPara.init_objects) {
            update_kDNN_object_multi(1, object_node, k_max_for_toain, k, mems.dist, mems.visited, mems.q,
                                     hier_local_knn_arr_single);
        }
        else{
            update_kDNN_object_multi(0, object_node, k_max_for_toain, k, mems.dist, mems.visited, mems.q,
                                     hier_local_knn_arr_single);
        }
//        cout<<"time: "<<clock()-start<<endl;
    }

}


#endif //TOAIN_FUNCTIONTEMPLATE_H
