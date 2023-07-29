#ifndef __METHOD_H__
#define __METHOD_H__

#include "query.h"

struct Method {
    virtual void index();
    virtual void save_index(const char *filename);
    virtual void load_index(const char *filename);

    virtual double query(const DistanceQuery &q) = 0;
};

#endif