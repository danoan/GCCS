#ifndef SEGBYCUT_TRIVIALCHECKER_H
#define SEGBYCUT_TRIVIALCHECKER_H

#include "MarkedMapChecker/MarkedMapCheckerInterface.h"

template<typename CheckableElement>
class TrivialChecker:public MarkedMapCheckerInterface<CheckableElement>
{
public:
    bool operator()(CheckableElement& sp) const{ return true;}

    void mark(CheckableElement& sp){}
    void unmark(CheckableElement& sp) {}
};

#endif //SEGBYCUT_TRIVIALCHECKER_H
