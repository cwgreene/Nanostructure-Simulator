#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define sq(x) ((x)*(x))

typedef double vector3[3];
typedef double vector3p[4]; //Vector plus data

typedef struct kdtree3
{
	vector3 location;
	int id; //associated data, is id in application
	int axis;
	struct kdtree3 *leftChild;
	struct kdtree3 *rightChild;
}kdtree3;

//Macro that generates comaprison function
#define COMPARE(direction,index) \
int compare3##direction (const vector3p a, const vector3p b)\
{\
	if(a[index] < b[index])	\
		return -1;\
	if(a[index] > b[index])\
		return 1;\
	return 0;\
}


COMPARE(x,0);
COMPARE(y,1);
COMPARE(z,2);
/*
int comparex(const vector3p a, const vector3p b)
{
	if(a[0] < b[0])	
		return -1;
	if(a[0] > b[0])
		return 1;
	return 0;
}

int comparey(const vector3p a, const vector3p b)
{
	if(a[1] < b[1])	
		return -1;
	if(a[1] > b[1])
		return 1;
	return 0;
}

int comparez(const vector3p a, const vector3p b)
{
	if(a[2] < b[2])	
		return -1;
	if(a[2] > b[2])
		return 1;
	return 0;
}*/



typedef int (*comparisonFunc3)(const vector3p a,const vector3p b);

comparisonFunc3 compare3[3] = {compare3x,compare3y,compare3z};

kdtree3 *kdtree3_cons()
{
	kdtree3 *newkd = (kdtree3 *)malloc(sizeof(kdtree3));
	newkd->leftChild = NULL;
	newkd->rightChild = NULL;
	return newkd;
}

kdtree3 *new_kdtree3(vector3p *array,int length,int depth)
{
	int axis = depth % 3;
	int median;
	//printf("Length:%d, Depth:%d, length/2:%d\n",length,depth,length/2);
	if(length <=0)
		return NULL;
	qsort ((void *)array, 
		length, sizeof(vector3p), 
		(int(*)(const void*,const void*))compare3[axis]);
	median = length/2;

	kdtree3 *kdTree = kdtree3_cons();
	kdTree->id = (int)array[median][3];//Won't work, need to attach id
//	printf("%d,%lf,%lf\n",kdTree->id,
//				array[median][0],
//				array[median][1]);
	kdTree->location[0] = array[median][0];
	kdTree->location[1] = array[median][1];
	kdTree->location[2] = array[median][2];
	kdTree->axis = axis;
	kdTree->leftChild = new_kdtree3(array,
				      median, depth+1);
	kdTree->rightChild = new_kdtree3(array+median+1,
				       (length-median)-1,depth+1);

	return kdTree;
}

//double dot3_c(double a1,double a2,double a3,double b1,double b2,double b3)
//{
//	return a1*b1+b2*a2+a3*b3;
//}

double dist3(vector3 *a,vector3 *b)
{
	//printf("dist2:\na:%lf,%lf\nb:%lf,%lf\n",(*a)[0],(*a)[1],
	//		(*b)[0],(*b)[1]);
	double a1 = (*a)[0];
	double a2 = (*a)[1];
	double a3 = (*a)[2];
	double b1 = (*b)[0];
	double b2 = (*b)[1];
	double b3 = (*b)[2];
	double result = ( sq(a1-b1)+sq(a2-b2)+sq(a3-b3) );
	//printf("dist_is:%lf\n\n",result);
	return result;
}

vector3 *find_point3_r(vector3 *point,kdtree3 *node,
			vector3 *best,double *bdist)
{
	int axis = node->axis;
	double cur_dist= dist3(point,&(node->location));
	kdtree3 *other_branch = NULL;
	//printf("point:%lf,%lf\n",(*point)[0],(*point)[1]);
	//printf("Visiting:%lf,%lf best:%lf\n",
		//node->location[0],node->location[1],
		//*bdist);

	//if cur_dist is smaller than the current best
	if(*bdist > cur_dist) 
	{
		(*best)[0] = node->location[0];
		(*best)[1] = node->location[1];
		(*best)[2] = node->location[2];
		//printf("new best: %lf,%lf\n",(*best)[0],(*best)[1]);
		
		*bdist = cur_dist;
	}
	
	if((*point)[axis] < node->location[axis]){
		if(node->leftChild != NULL)
			best = find_point3_r(point,node->leftChild,best,bdist);
		other_branch= node->rightChild;
	}
	else{
		if(node->rightChild != NULL)
			best = find_point3_r(point,node->rightChild,best,bdist);
		other_branch= node->leftChild;
	}

	//How close are we to other side?
	double split_dist = (*point)[axis] - node->location[axis];
	split_dist = split_dist*split_dist; //bdist is squared
	if(split_dist<*bdist) //could be closer on other side
	{	
		if(other_branch != NULL)
			best = find_point3_r(point,other_branch,best,bdist);
	}
	return best;
}
static int static_count = 0;
vector3 *find_point3_r_id(vector3 *point,kdtree3 *node,
			vector3 *best,double *bdist,int *id)
{
	int axis = node->axis;
	double cur_dist= dist3(point,&(node->location));
	kdtree3 *other_branch = NULL;
	//printf("point:%lf,%lf\n",(*point)[0],(*point)[1]);
	//printf("Visiting:%lf,%lf best:%lf\n",
	//	node->location[0],node->location[1],
	//	*bdist);
	
	//if( (static_count++) % (1000*1000) == 0)
	//	printf("called: %d \n",static_count);

	//if cur_dist is smaller than the current best
	if(*bdist > cur_dist) 
	{
		(*best)[0] = node->location[0];
		(*best)[1] = node->location[1];
		(*best)[2] = node->location[2];
		//printf("new best: %lf,%lf id:%id\n",(*best)[0],(*best)[1],
				//node->id);
		*id = node->id;
		*bdist = cur_dist;
	}
	
	if((*point)[axis] < node->location[axis]){
		if(node->leftChild != NULL)
			best = find_point3_r_id(point,node->leftChild,best,bdist,id);
		other_branch= node->rightChild;
	}
	else{
		if(node->rightChild != NULL)
			best = find_point3_r_id(point,node->rightChild,best,bdist,id);
		other_branch= node->leftChild;
	}

	//How close are we to other side?
	double split_dist = (*point)[axis] - node->location[axis];
	split_dist = split_dist*split_dist; //bdist is squared
	if(split_dist<*bdist) //could be closer on other side
	{	
		if(other_branch != NULL)
			best = find_point3_r_id(point,other_branch,
						best,bdist,id);
	}
	return best;
}

int kdtree_find_point3_id(kdtree3 *tree,vector3 *point)
{
	vector3 *best = malloc(sizeof(vector3));
	(*best)[0] = 100000.0;
	(*best)[1] = 100000.0;
	(*best)[2] = 100000.0;
	//long unsigned int start = clock();
	int id= 0;
	double bdist = dist3(point,best);
	//printf("bdist: %lf",bdist);
	find_point3_r_id(point,tree,best,&bdist,&id);
	//printf("return id: %d\n",id);
	//printf("elapsed: %lf\n",((double) ((clock()-start))/CLOCKS_PER_SEC));
	free(best);
	return id;
}

//Bad function! Will have memory links
vector3 *kdtree_find_point3(kdtree3 *tree,vector3 *point)
{
	vector3 *best = malloc(sizeof(vector3));
	(*best)[0] = 100000.0;
	(*best)[1] = 100000.0;
	(*best)[2] = 100000.0;
	printf(
	"You really should not use kdtree_find_point3, has Memory leaks\n");
	//printf("starting point:%lf,%lf\n",(*point)[0],(*point)[1]);
	double bdist = dist3(point,best);
	best = find_point3_r(point,tree,best,&bdist);
	return best;
}

void print_kdtree3(kdtree3 *tree)
{
	if(tree == NULL)
		return;
	printf("%lf,%lf,%lf\n",tree->location[0],
				tree->location[1],
				tree->location[2]);
	print_kdtree3(tree->leftChild);
	print_kdtree3(tree->rightChild);
}

void call_this3()
{
	printf("Hello Python!\n");
}

int do_nothing3(int x)
{
	return x;
}
