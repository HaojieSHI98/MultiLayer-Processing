/*
 * Point.h
 *
 *  Created on: 2017年2月14日
 *      Author: lostrong
 */

#ifndef SRC_POINT_H_
#define SRC_POINT_H_

struct Point
{
	int x;
	int y;
	int flag;
	Point() :x(0), y(0), flag(0){}
	Point(int a, int b, int c) :x(a), y(b), flag(c){}
	static bool compare_x(const Point&a, const Point&b){ return a.x < b.x; }
	static bool compare_y(const Point&a, const Point&b){ return a.y < b.y; }
};



#endif /* SRC_POINT_H_ */
