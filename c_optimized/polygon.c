typedef double vector2[2];

#define dot(x,y) (x[0]*y[0]+x[1]*y[1])
#define COPY(vec1,vec2) vec1[0] = vec2[0]; vec1[1] = vec2[1];
#define MINUS(x,y,dest) dest[0] = x[0] - y[0]; dest[1] = x[1]-y[1];
#define SCALE(x,scale,dest) dest[0] = x[0]/scale; dest[1] = x[1]/scale;

typedef struct
{
	vector2 *normals;
	vector2 *midpoints;
	vector2 *vertices;
}polygon;

polygon *polygon_from_points(vector2 *points,int length)
{
	int i = 0;
	vector2 result;
	polygon *poly = malloc(sizeof(polygon));
	poly->vertices = malloc(sizeof(vector2)*length);
	poly->normals  = malloc(sizeof(vector2)*length);
	poly->midpoints = malloc(sizeof(vector2)*length);

	for(i=0;i<length;i++)
	{
		COPY(poly->vertices[i],points[i]);
		//midpoints
		poly->midpoints[i][0] = points[i][0]-points[(i+1)%length][1];
		poly->midpoints[i][1] = points[i][0]-points[(i+1)%length][1];
		//normals
	}
}

int inside(vector2 *normals, vector2 *midpoint,
	   int length,vector2 point)
{
	int i;
	vector2 result;
	for(i = 0; i < length;i++)
	{
		//midpoint[i] - point
		MINUS(midpoint[i],point,result);
		if(dot(normals[i],result) < 0)
			return false; //outside
	}
	return true;//inside
}
