#ifndef DEBUG_H
#define DEBUG_H
#include <iostream>
namespace debug{
class DebugOut : public std::ostream
{
};
DebugOut dbg;
}
#define TRACK //debug tracking code
#ifdef TRACK
static int grabbed = -1;
#endif 
#endif 
