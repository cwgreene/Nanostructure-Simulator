#ifndef RAND_HPP
#define RAND_HPP
#include <stdlib.h>

//randint(int min, int max)
//inclusively generate between min and max
int randint(int a, int b)
{
	int range = (b-a)+1;
	return ((rand()%range)+a);
}
#endif
