//
// Created by lostrong on 17/7/25.
//

#ifndef SOB_UCARINFO_H
#define SOB_UCARINFO_H
class UCarInfo{
public:
    long id;
    long timestamp;
    double x;
    double y;
    double speed;
    int angle;
    UCarInfo(long id, long timestamp, double x, double y, double speed, int angle){
        this->id=id;
        this->timestamp=timestamp;
        this->x=x;
        this->y=y;
        this->speed=speed;
        this->angle=angle;
    }
};
#endif //SOB_UCARINFO_H
