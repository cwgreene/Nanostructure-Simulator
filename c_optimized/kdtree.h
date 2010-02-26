#ifndef KDTREE_H
#define KDTREE_H
typedef double vector2[2];
typedef double vector2p[2];

typedef struct kdtree
{
	vector2 location;
	int axis;
	struct kdtree *leftChild;
	struct kdtree *rightChild;
}kdtree;

kdtree *new_kdtree(vector2 *array,int length,int depth);
void print_kdtree(kdtree *tree);
vector2 *kdtree_find_point(kdtree *tree,vector2 *point,vector2 *best);
int kdtree_find_point_id(kdtree *tree,vector2 *point);
#endif
