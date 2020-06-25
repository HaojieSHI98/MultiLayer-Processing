//
// Created by lostrong on 17/8/4.
//

#ifndef SOB_CARPERIODLOCATION_H
#define SOB_CARPERIODLOCATION_H

class CarPeriodLocation{
public:
    int car_id;
    int period_num;
    int nearest_node_id;
    CarPeriodLocation(int car_id, int period_num, int nearest_node_id){
        this->car_id=car_id;
        this->period_num=period_num;
        this->nearest_node_id=nearest_node_id;
    }

};
#endif //SOB_CARPERIODLOCATION_H
