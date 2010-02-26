#include "kdtree.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
	vector2p bob[3];
	bob[0][0] = 0.;
	bob[0][1] = 1.;
	bob[0][2] = 0;
	bob[1][0] = 2.;
	bob[1][1] = .3;
	bob[1][2] = 2;
	bob[2][0] = -.2;
	bob[2][1] = 1.;
	bob[2][2] = 3;

	printf("hi!\n");	
	kdtree *joe= new_kdtree(bob,3,0);	
	print_kdtree(joe);

	vector2 closest;
	printf("bob[0]: %lf,%lf\n",bob[0][0],bob[0][1]);
	kdtree_find_point(joe,&(bob[0]),&closest);
	printf("Closest: %lf,%lf\n",closest[0],closest[1]);
	return 0;
}
