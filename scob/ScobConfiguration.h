//
// Created by lostrong on 17/7/13.
//

#ifndef SOB_SCOBCONFIGURATION_H
#define SOB_SCOBCONFIGURATION_H
#include <vector>

using namespace std;
class ScobConfiguration{
public:
    int cut_level;
    int tradeoff_level;
    int mode;
    ScobConfiguration(int tradeoff_level, int cut_level, int mode){
        this->cut_level=cut_level;
        this->tradeoff_level=tradeoff_level;
        this->mode=mode;
    }
    static vector<ScobConfiguration> generate_all_configurations(int level_num) {
        vector<ScobConfiguration> configurations;
        int tradeoff_level = 0;
        for (int cut_level = 0; cut_level <= level_num; cut_level++) {
            configurations.push_back(*(new ScobConfiguration(tradeoff_level, cut_level, QUERY_FAVOR_MODE)));
        }
        tradeoff_level = level_num;
        for (int cut_level = level_num; cut_level >= 0; cut_level--) {
            configurations.push_back(*(new ScobConfiguration(tradeoff_level, cut_level, UPDATE_FAVOR_MODE)));
        }
        return configurations;
    }

    static vector<ScobConfiguration> generate_all_update_favor_configurations(int level_num) {
        vector<ScobConfiguration> configurations;
        int tradeoff_level = level_num;
        for (int cut_level = level_num; cut_level >= 0; cut_level--) {
            configurations.push_back(*(new ScobConfiguration(tradeoff_level, cut_level, UPDATE_FAVOR_MODE)));
        }
        return configurations;
    }

};
#endif //SOB_SCOBCONFIGURATION_H
