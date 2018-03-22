#ifndef SEGBYCUT_LAZYCOMBINATIONS_H
#define SEGBYCUT_LAZYCOMBINATIONS_H

#include <stack>
#include "PropertyChecker/MarkedMapChecker/MarkedMapCheckerInterface.h"
#include "PropertyChecker/TrivialChecker.h"

template<typename Container, int HowMany>
class LazyCombinations
{
public:
    enum OperationType{Offspring,Sibling,SiblingAndUnmark,Done,Invalid};

    class Indications
    {
    public:
        Indications(){};

        Indications(int level,
                    int index,
                    int finalIndex,
                    OperationType operation):level(level),
                                             index(index),
                                             finalIndex(finalIndex),
                                             operation(operation)
        {
            assertSoundness();
        }

        void assertSoundness()
        {
            if(level==HowMany) this->operation = Done;
            else if(index>finalIndex) this->operation = Invalid;
        }

    public:
        int level;
        int index;

        int finalIndex;
        OperationType operation;
    };

    typedef typename Container::value_type ContainerValueType;

    BOOST_CONCEPT_ASSERT( (CheckableElementConcept<ContainerValueType>) );

public:
    LazyCombinations(Container& container):_container(container)
    {
        _stack.push( Indications(0,0,container.size()-HowMany,Offspring) );
    }

    ~LazyCombinations(){ /*delete _checker;*/ }

    void addConsistencyChecker(MarkedMapCheckerInterface<ContainerValueType>* checker)
    {
        _checkersList.push_back( checker );
    }

    bool next(ContainerValueType* nextElement)
    {
        Indications i;
        do
        {
            i = _stack.top();
            _stack.pop();


            switch(i.operation)
            {
                case Offspring:
                    if( checkElement( _container[i.index] ) )  //This seed pair intersects with a previous seed pair
                    {
                        i.operation = SiblingAndUnmark;
                        _stack.push(i);

                        markElement(_container[i.index]);

                        _stack.push(Indications(i.level + 1,
                                                i.index + 1,
                                                i.finalIndex + 1,
                                                Offspring));

                        _current[i.level] = i.index;
                    }
                    else
                    {
                        i.operation = Sibling;
                        _stack.push(i);
                    }


                    break;
                case SiblingAndUnmark:
                    unmarkElement(_container[i.index]);
                case Sibling:
                    i.index += 1;
                    i.operation = Offspring;
                    i.assertSoundness();

                    _stack.push(i);
                    break;
            }


        }while(!_stack.empty() && i.operation!=Done);

        if(!_stack.empty())
        {
            //Done, return next combination
            for(int j=0;j<HowMany;++j)
            {
                nextElement[j] = _container[ _current[j] ].data();
            }

            return true;
        }
        else
        {
            return false;
        }
    }

private:
    bool checkElement(ContainerValueType& v)
    {
        for(auto it=_checkersList.begin();it!=_checkersList.end();++it)
        {
            if( !(*it)->operator()(v) ) return false;
        }
        return true;
    }

    void markElement(ContainerValueType& v)
    {
        for(auto it=_checkersList.begin();it!=_checkersList.end();++it)
        {
            (*it)->mark(v);
        }
    }

    void unmarkElement(ContainerValueType& v)
    {
        for(auto it=_checkersList.begin();it!=_checkersList.end();++it)
        {
            (*it)->unmark(v);
        }
    }

public:
    Container& _container;
    std::stack<Indications> _stack;

    std::vector< MarkedMapCheckerInterface< ContainerValueType >* > _checkersList;

    int _current[HowMany];

};






#endif //SEGBYCUT_LAZYCOMBINATIONS_H
