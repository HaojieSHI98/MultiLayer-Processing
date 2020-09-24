/*
 * GlobalVariables.h
 *
 *  Created on: 2016年7月6日
 *      Author: lostrong
 */

#ifndef SRC_ROADGLOBALVARIABLES_H_
#define SRC_ROADGLOBALVARIABLES_H_

#include<vector>
#include<string>
#include<unordered_map>
#include<limits.h>
//#include<omp.h>
#include <thread>
#include "../graph/CellNode.h"
#include "../object/CarPos.h"
#include "../object/CarPeriodLocation.h"
#include "TestParameter.h"

using namespace std;

double EARTH_RADIUS = 6378.137; //地球半径
double PI = 3.1415926535;


//3000000 for NY, FLA
//7000000 for W
// MAGIC_NUM 5000000 for NY, FLA
// MAGIC_NUM 10000000 for W
#define MAGIC_NUM 2000000
#define NODE_VERSION_RAN 1
#define NODE_VERSION_REAL 2
#define EDGE_VERSION_RAN 3
#define EDGE_VERSION_REAL 4
#define MICROSEC_PER_SEC 1000000
#define MILLISEC_PER_SEC 1000
#define MAX_LEVELS 12
#define NUM_CELL_IN_NEIGHBOR 25
#define QUERY_FAVOR_MODE 2
#define UPDATE_FAVOR_MODE 3
#define QUEUE_POLICY_QUERY_FIRST 1
#define QUEUE_POLICY_FIFO 0
#define DEFAULT_K 20
#define DEFAULT_RATIO 0.001
#define DATA_RANDOM 0
#define DATA_BJ_REAL 1
#define DATA_NW_REAL 2
#define SHUFFLE_NUM 1000000
#define BJ_OBJECTS 2770
#define NW_OBJECTS 13000
#define ZIPF_OBJECT 0
#define ZIPF_QUERY 0
#define NUM_OBV_T 1
//#define DISPLAY 0
//new

#define NORMAL_MODE 0
#define EVALUATION_START_MODE 1
#define EVALUATION_END_MODE 2
#define RESET_MODE 3
#define UPDATE_SET 0
#define QUERY_SET 1
#define MIN_UQ_DIFF 10
#define EVA_TIME 5

#define EXP_SIZE 1000
#define STAR_NUM 20
#define X_STAR_MODE 0
#define FILTER_NUM 5
#define MAX_NORECORD_TIME 2
int No_Record_Flag=0;
long Start_No_Record_Time = 0;
int Update_Query_Threshold = 1;
typedef struct{
    std::mutex ta_mutex[2];
    std::mutex tq_mutex[2];
    std::mutex tu_mutex[2];
    double query_rate=0;
    double update_rate=0;
    vector<int> task_list;
    vector<double> time_list;
    double ratio_x[2]={0};
    vector<double> ts[2];
    vector<double> ta[2];
    vector<double> tq_ex[2];
    vector<double> tu_ex[2];
    double tq[2]={0};
    double tu[2]={0};
    double Vq[2]={0};
    double Vu[2]={0};
    double update_query_ratio=0;
    double last_update_query_ratio=0;
    vector<double> update_query_ratio_set;
}Observer;
vector<int> x_stars;
Observer observer;
int DISPLAY = 1;
//
int NEIGHBOR_SEARCH=0;
double ALPHA=2.0;
int COVERAGE_INDEX=0;
int BETWEEN_INDEX=0;
int ADA_COVERAGE_INDEX=0;
int ADA_BETWEEN_INDEX=0;
int verifyResult = 0;

int *x;
int *y;
vector<int> **cells;
vector<int> **circle_cells;
int *cell_loc;
int *circle_cell_loc[25];
int eNext[2 * MAGIC_NUM + 2] = {0}, eDis[2 * MAGIC_NUM + 1] = {0}, eNode[2 * MAGIC_NUM + 1] = {0};
unordered_map<int, int> valid_node_hash;
int *vStart;
int eNum = 1;
int test_n;
int test_m;
int cell_num;

int *regional_landmarks_all;
double largest_distance;
unordered_map<string, string> data_map;
unordered_map<string, string> network_map;
string data_space;
int *node_to_cell;

//Parameters to control the algorithm
int max_level = 20;
int tradeoff_level = 0;
int cut_level = 9;
int mode = 2;
int level_num = 0;

//The number of nodes within a leaf
int leaf_num = 50;
//int leaf_num = 100;

//Global indicator for data conflicts
int has_read_road_network = 0;
int has_read_scob_index = 0;
int has_read_ucar_data=0;
int max_period_num=0;
int memory_deleted=0;
unordered_map<int, vector<KNode> > shortcuts;
vector<KNode> *hier_local_knn_arr;
vector<KNode> **level_scs_arr;
vector<KNode> **shortcuts_arr_local;
vector<KNode> **shortcuts_inc_arr_local;
vector<KNode> shortcuts_arr[MAGIC_NUM];
vector<KNode> *rev_shortcuts_arr;
vector<KNode> *rev_shortcuts_level_arr;
vector<KNode> shortcuts_inc_arr[MAGIC_NUM];
vector<int> car_localarr_map[30000];
int *node_order = NULL;
int *inlev_order = NULL;
int cell_nums[50];
int level_gap[17] = {2, 5, 11, 23, 47, 95, 191, 383, 767, 1535, 3071, 6143, 12287, 24575, 49151, 98303, 196607};
int *cell_loc_x;
int *cell_loc_y;
vector<std::pair<int, int> > node_ranks;
int x_min = INT_MAX, y_min = INT_MAX, x_max = INT_MIN, y_max = 0, x_len = 0, y_len = 0;



TestParameter test_parameters;
InputParameter input_parameters;
//string network_name="BJ-old";
string network_name="NY";


unordered_map<int, CarPos> carposes_map;
vector<CarPos> carposes;
vector<CarPos> carposes_random;
vector<CarPos> allcarposes;
//node version, random
vector<int> objects_vector;
unordered_map<int, vector<PointCar> > pointcar_hash;
unordered_map<int, vector<PointCar> > edge_car_hash;

// ucar information
vector<vector<CarPeriodLocation> > car_period_location_list;
vector<int> indexes;
vector<int> valid_cars;



// just a static block, for initialize certain global variables;
struct StaticBlock {
    StaticBlock() {
        data_map["bj1_node.csv"] = "0";
//        data_map["USA-road-d.BJ.co"] = "0";
        data_map["USA-road-d.NY.co"] = "1";
        data_map["USA-road-d.FLA.co"] = "2";
        data_map["USA-road-d.LKS.co"] = "3";
        data_map["USA-road-d.NE.co"] = "4";
        data_map["USA-road-d.W.co"] = "5";
        data_map["USA-road-d.USA.co"] = "6";
        data_map["USA-road-d.NW.co"] = "7";
        data_map["USA-road-d.BJ10000.co"] = "8";
        data_map["USA-road-d.BJ20000.co"] = "9";
        data_map["USA-road-d.BJ40000.co"] = "10";
        data_map["USA-road-d.BJ80000.co"] = "11";
        data_map["USA-road-d.BJ160000.co"] = "12";
        data_map["USA-road-d.BJ320000.co"] = "13";
        data_map["USA-road-d.BJ640000.co"] = "14";
        data_map["USA-road-d.BJ1280000.co"] = "15";
        network_map["bj1_node.csv"] = "BJ-old";
            network_map["USA-road-d.BJ.co"] = "BJ";
        network_map["USA-road-d.NY.co"] = "NY";
        network_map["USA-road-d.FLA.co"] = "FLA";
        network_map["USA-road-d.LKS.co"] = "LKS";
        network_map["USA-road-d.NE.co"] = "NE";
        network_map["USA-road-d.W.co"] = "W";
        network_map["USA-road-d.USA.co"] = "USA";
        network_map["USA-road-d.BJ10000.co"] = "BJ10000";
        network_map["USA-road-d.BJ20000.co"] = "BJ20000";
        network_map["USA-road-d.BJ40000.co"] = "BJ40000";
        network_map["USA-road-d.BJ80000.co"] = "BJ80000";
        network_map["USA-road-d.BJ160000.co"] = "BJ160000";
        network_map["USA-road-d.BJ320000.co"] = "BJ320000";
        network_map["USA-road-d.BJ640000.co"] = "BJ640000";
        network_map["USA-road-d.BJ1280000.co"] = "BJ1280000";
    }

    ~StaticBlock(){

    }
};


static StaticBlock staticBlock;

#endif /* SRC_ROADGLOBALVARIABLES_H_ */
