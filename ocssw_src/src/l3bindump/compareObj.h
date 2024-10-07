/* 
 * File:   compareObj.h
 * Author: dshea
 *
 * Created on July 16, 2015, 2:33 PM
 */

#ifndef COMPAREOBJ_H
#define COMPAREOBJ_H

class CompareObj {
public:
    virtual bool isInside(double lat, double lon) = 0;
};

class CompareObjLonOnly : public CompareObj {
public:
    CompareObjLonOnly(double east, double west);

    virtual bool isInside(double lat, double lon);

private:
    double east;
    double west;
};

class CompareObjCircle : public CompareObj {
public:
    CompareObjCircle(double lat, double lon, double radius);

    virtual bool isInside(double lat, double lon);

private:
    double lat0;
    double lon0;
    double radius;
};

#endif /* COMPAREOBJ_H */

