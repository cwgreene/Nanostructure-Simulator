#ifndef KDTREE3_H
#define KDTREE3_H
typedef double vector3[3];
typedef double vector3p[3];

typedef struct kdtree3
{
	vector3 location;
	int axis;
	struct kdtree *leftChild;
	struct kdtree *rightChild;
}kdtree3;

kdtree3 *new_kdtree3(vector3 *array,int length,int depth);
void print_kdtree3(kdtree3 *tree);
vector3 *kdtree_find_point3(kdtree3 *tree,vector3 *point,vector3 *best);
int kdtree_find_point3_id(kdtree3 *tree,vector3 *point);
#endif
