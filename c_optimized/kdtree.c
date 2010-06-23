#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef double vector2[2];
typedef double vector2p[3];

typedef struct kdtree
{
	vector2 location;
	int id;
	int axis;
	struct kdtree *leftChild;
	struct kdtree *rightChild;
}kdtree;


int comparex(const vector2p a, const vector2p b)
{
	if(a[0] < b[0])	
		return -1;
	if(a[0] > b[0])
		return 1;
	return 0;
}

int comparey(const vector2p a, const vector2p b)
{
	if(a[1] < b[1])	
		return -1;
	if(a[1] > b[1])
		return 1;
	return 0;
}

typedef int (*comparisonFunc)(const vector2p a,const vector2p b);

comparisonFunc compare[2] = {comparex,comparey};

kdtree *kdtree_cons()
{
	kdtree *newkd = (kdtree *)malloc(sizeof(kdtree));
	newkd->leftChild = NULL;
	newkd->rightChild = NULL;
	return newkd;
}

kdtree *new_kdtree(vector2p *array,int length,int depth)
{
	int axis = depth % 2;
	int median;
	//printf("Length:%d, Depth:%d, length/2:%d\n",length,depth,length/2);
	if(length <=0)
		return NULL;
	qsort ((void *)array, 
		length, sizeof(vector2p), 
		(int(*)(const void*,const void*))compare[axis]);
	median = length/2;

	kdtree *kdTree = kdtree_cons();
	kdTree->id = (int)array[median][2];//Won't work, need to attach id
//	printf("%d,%lf,%lf\n",kdTree->id,
//				array[median][0],
//				array[median][1]);
	kdTree->location[0] = array[median][0];
	kdTree->location[1] = array[median][1];
	kdTree->axis = axis;
	kdTree->leftChild = new_kdtree(array,
				      median, depth+1);
	kdTree->rightChild = new_kdtree(array+median+1,
				       (length-median)-1,depth+1);

	return kdTree;
}

double dot_c(double a1,double a2,double b1,double b2)
{
	return a1*b1+b2*a2;
}

double dist2(vector2 *a,vector2 *b)
{
	//printf("dist2:\na:%lf,%lf\nb:%lf,%lf\n",(*a)[0],(*a)[1],
	//		(*b)[0],(*b)[1]);
	double a1 = (*a)[0];
	double a2 = (*a)[1];
	double b1 = (*b)[0];
	double b2 = (*b)[1];
	double result = ( (a1-b1)*(a1-b1)+
	   	 (a2-b2)*(a2-b2));
	//printf("dist_is:%lf\n\n",result);
	return result;
}

vector2 *find_point_r(vector2 *point,kdtree *node,
			vector2 *best,double *bdist)
{
	int axis = node->axis;
	double cur_dist= dist2(point,&(node->location));
	kdtree *other_branch = NULL;
	//printf("point:%lf,%lf\n",(*point)[0],(*point)[1]);
	//printf("Visiting:%lf,%lf best:%lf\n",
		//node->location[0],node->location[1],
		//*bdist);

	//if cur_dist is smaller than the current best
	if(*bdist > cur_dist) 
	{
		(*best)[0] = node->location[0];
		(*best)[1] = node->location[1];
		//printf("new best: %lf,%lf\n",(*best)[0],(*best)[1]);
		
		*bdist = cur_dist;
	}
	
	if((*point)[axis] < node->location[axis]){
		if(node->leftChild != NULL)
			best = find_point_r(point,node->leftChild,best,bdist);
		other_branch= node->rightChild;
	}
	else{
		if(node->rightChild != NULL)
			best = find_point_r(point,node->rightChild,best,bdist);
		other_branch= node->leftChild;
	}

	//How close are we to other side?
	double split_dist = (*point)[axis] - node->location[axis];
	split_dist = split_dist*split_dist; //bdist is squared
	if(split_dist<*bdist) //could be closer on other side
	{	
		if(other_branch != NULL)
			best = find_point_r(point,other_branch,best,bdist);
	}
	return best;
}
//static int static_count = 0;
vector2 *find_point_r_id(vector2 *point,kdtree *node,
			vector2 *best,double *bdist,int *id)
{
	int axis = node->axis;
	double cur_dist= dist2(point,&(node->location));
	kdtree *other_branch = NULL;
	//printf("point:%lf,%lf\n",(*point)[0],(*point)[1]);
	//printf("Visiting:%lf,%lf best:%lf\n",
		//node->location[0],node->location[1],
		//*bdist);
	
	//if( (static_count++) % (1000*1000) == 0)
		//printf("called: %d \n",static_count);

	//if cur_dist is smaller than the current best
	if(*bdist > cur_dist) 
	{
		(*best)[0] = node->location[0];
		(*best)[1] = node->location[1];
		//printf("new best: %lf,%lf\n",(*best)[0],(*best)[1]);
		*id = node->id;
		*bdist = cur_dist;
	}
	
	if((*point)[axis] < node->location[axis]){
		if(node->leftChild != NULL)
			best = find_point_r_id(point,node->leftChild,best,bdist,id);
		other_branch= node->rightChild;
	}
	else{
		if(node->rightChild != NULL)
			best = find_point_r_id(point,node->rightChild,best,bdist,id);
		other_branch= node->leftChild;
	}

	//How close are we to other side?
	double split_dist = (*point)[axis] - node->location[axis];
	split_dist = split_dist*split_dist; //bdist is squared
	if(split_dist<*bdist) //could be closer on other side
	{	
		if(other_branch != NULL)
			best = find_point_r_id(point,other_branch,best,bdist,id);
	}
	return best;
}

int kdtree_find_point_id(kdtree *tree,vector2 *point)
{
	vector2 *best = malloc(sizeof(vector2));
	(*best)[0] = 100000.0;
	(*best)[1] = 100000.0;
	//long unsigned int start = clock();
	int id= 0;
	double bdist = dist2(point,best);
	find_point_r_id(point,tree,best,&bdist,&id);
	//printf("elapsed: %lf\n",((double) ((clock()-start))/CLOCKS_PER_SEC));	
	free(best);
	return id;
}


/*Very bad function! Will have memory leaks.*/
vector2 *kdtree_find_point(kdtree *tree,vector2 *point)
{
	vector2 *best = malloc(sizeof(vector2));
	(*best)[0] = 100000.0;
	(*best)[1] = 100000.0;
	//printf("starting point:%lf,%lf\n",(*point)[0],(*point)[1]);
	double bdist = dist2(point,best);
	best = find_point_r(point,tree,best,&bdist);
	return best;
}

void print_kdtree(kdtree *tree)
{
	if(tree == NULL)
		return;
	printf("%lf,%lf\n",tree->location[0],tree->location[1]);
	print_kdtree(tree->leftChild);
	print_kdtree(tree->rightChild);
}

#define max(x,y) (((x) > (y) ? (x) : (y)))
int longest_chain(kdtree *tree,int depth,int maxdepth)
{
	if(tree->leftChild ==NULL && tree->rightChild == NULL)
	{
		return max(depth,maxdepth);
	}
	if(tree->rightChild != NULL)
	{
		maxdepth = longest_chain(tree->rightChild,depth+1,maxdepth);
	}
	if(tree->leftChild != NULL)
	{
		maxdepth = longest_chain(tree->leftChild,depth+1,maxdepth);
	}
	return maxdepth;
}

void call_this()
{
	printf("Hello Python!\n");
}

int do_nothing(int x)
{
	return x;
}
