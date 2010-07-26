#ifndef DEBUG_H
#define DEBUG_H
#include <iostream>
namespace debug{
class DebugOut : public std::ostream
{
};
DebugOut dbg;
}
#endif 
