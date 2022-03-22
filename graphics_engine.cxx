//My movie can be found at:
//https://youtu.be/rlv7rD3e_bU

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <cmath>
#include <string>
#include <algorithm>

#define NORMALS
#define WIDTH 1000
#define HEIGHT 1000

using namespace std;

void Normalize(double *point)
{	
	double magnitude = sqrt(point[0]*point[0]+point[1]*point[1]+point[2]*point[2]);
	for(int i=0; i<3; ++i)
	{
		point[i] = point[i]/magnitude;
	}
}

double *Normalize(double x, double y, double z)
{
	double *vec = new double[3];
	vec[0] = x;
	vec[1] = y;
	vec[2] = z;

	Normalize(vec);

	return vec;
}

class Point
{
	public:
		double colors[3];
		double x;
		double y;
		double z;
		double w = 1;
		double shading;
		double normal[3];
};

double ceil__mod(double f)
{
    return ceil(f-0.00001);
}

double floor__mod(double f)
{
    return floor(f+0.00001);
}


vtkImageData *
NewImage(int height, int width)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}


void CrossProduct(double *left, double *right, double *result)
{
		
	result[0] = left[1]*right[2] - left[2]*right[1];
	result[1] = left[2]*right[0] - left[0]*right[2];
	result[2] = left[0]*right[1] - left[1]*right[0];

}

double DotProduct(double *left, double *right)
{
	double result=0;
	for(int i=0; i<3; ++i)
	{
		result+= left[i]*right[i];
	}
	
	return result;
}

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const Point *ptIn, Point *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
    {
		  for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }
    }
	 return rv;
}


void
Matrix::TransformPoint(const Point *ptIn, Point *ptOut)
{
    ptOut->x = ptIn->x*A[0][0]
             + ptIn->y*A[1][0]
             + ptIn->z*A[2][0]
             + ptIn->w*A[3][0];
    
	 ptOut->y = ptIn->x*A[0][1]
             + ptIn->y*A[1][1]
             + ptIn->z*A[2][1]
             + ptIn->w*A[3][1];
    
	 ptOut->z = ptIn->x*A[0][2]
             + ptIn->y*A[1][2]
             + ptIn->z*A[2][2]
             + ptIn->w*A[3][2];
    
	 ptOut->w = ptIn->x*A[0][3]
             + ptIn->y*A[1][3]
             + ptIn->z*A[2][3]
             + ptIn->w*A[3][3];
}


class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          CameraTransform();
    Matrix          ViewTransform();
    Matrix          DeviceTransform();
};

Matrix Camera::CameraTransform()
{
//here we are making the camera frame

	//Calculate camera basis vectors
	double u[3], v[3], w[3], O[3];
	
	for(int i=0; i<3; ++i)
	{
		O[i] = position[i];
		w[i] = position[i] - focus[i] ;
		v[i] = up[i];
	}
	//calculate u
	CrossProduct(v, w, u);

	//use w and u to find v
	CrossProduct(w, u, v);
	
	Normalize(u);
	Normalize(v);
	Normalize(w);

//frame is made now we can cacluate transform values
	Matrix cameraTransform;
	cameraTransform.A[3][3] = 1;
	for(int i=0; i<3; ++i)
	{
		cameraTransform.A[i][3] = 0;
	}
	
	for(int i=0; i<3; ++i)
	{	
		int j=0;
		cameraTransform.A[i][j++] = u[i];
		cameraTransform.A[i][j++] = v[i];
		cameraTransform.A[i][j] = w[i];
	}
	
	double t[3]; //initialize t as (0,0,0) - O
	for(int i=0; i<3; ++i)
	{
		t[i] = 0 - O[i];
	}

	cameraTransform.A[3][0] = DotProduct(u,t);
	cameraTransform.A[3][1] = DotProduct(v,t);
	cameraTransform.A[3][2] = DotProduct(w,t);
	
	return cameraTransform;
}

Matrix Camera::ViewTransform()
{
	Matrix viewTransform;
	double cot = cos(angle/2)/sin(angle/2);
	for(int i=0; i<4; ++i)
	{
		for(int j=0; j<4; ++j)
		{
			viewTransform.A[i][j] = 0;
		}
	}
	viewTransform.A[0][0] = cot;
	viewTransform.A[1][1] = cot;
	viewTransform.A[2][2] = (far+near)/(far-near);
	viewTransform.A[2][3] = -1;
	viewTransform.A[3][2] = (2*far*near)/(far-near);

	return viewTransform;
}

Matrix Camera::DeviceTransform()
{
	Matrix deviceTransform;
	double n = WIDTH, m = HEIGHT;
	for(int i=0; i<4; ++i)
	{
		for(int j=0; j<4; ++j)
		{
			deviceTransform.A[i][j] = 0;
		}
	}
	
	deviceTransform.A[0][0] = n/2;
	deviceTransform.A[1][1] = m/2;
	deviceTransform.A[2][2] = 1;
	deviceTransform.A[3][0] = n/2;
	deviceTransform.A[3][1] = m/2;
	deviceTransform.A[3][3] = 1;
	
	return deviceTransform;
}



double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}


struct LightingParameters
{
    LightingParameters(void)
	 { 
	 	lightDir[0] = -0.6;
		lightDir[1] = 0;
	   lightDir[2] = -0.8;
      Ka = 0.3;
	   Kd = 0.7;
      Ks = 2.8;
      alpha = 50.5;
	 };
    double lightDir[3]; // The direction of the light source
    double Ka;          // The coefficient for ambient lighting
    double Kd;          // The coefficient for diffuse lighting
    double Ks;          // The coefficient for specular lighting
    double alpha;       // The exponent term for specular lighting
};


LightingParameters
GetLighting(Camera c)
{
    LightingParameters lp;
    lp.lightDir[0] = c.position[0]-c.focus[0];
    lp.lightDir[1] = c.position[1]-c.focus[1];
    lp.lightDir[2] = c.position[2]-c.focus[2];
    double mag = sqrt(lp.lightDir[0]*lp.lightDir[0]
                    + lp.lightDir[1]*lp.lightDir[1]
                    + lp.lightDir[2]*lp.lightDir[2]);
    if (mag > 0)
    {
        lp.lightDir[0] /= mag;
        lp.lightDir[1] /= mag;
        lp.lightDir[2] /= mag;
    }

    return lp;
}


class Triangle
{
  public:
      double         X[3];
      double         Y[3];
		double 			Z[3];
      double  colors[3][3];
		double normals[3][3];
		double shading[3];
  // would some methods for the triangle be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;
		double *zBuffer;
  // would some methods for accessing and setting pixels be helpful?
};

//This is the get Triangles method for 1D


std::vector<Triangle>
GetTriangles()
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("graphics_engine_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS

        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;

        }
    }

    return tris;
}



                          /*My functions and classes start here*/

void PrintPoints(Triangle triangle)
{	
	for(int i=0;i<3;++i)
	{
			cout<<i<<": ";
			cout<<"("<<triangle.X[i]<<","<<triangle.Y[i]<<")"<<endl;
	}
	cout<<endl;
}

//x is the point between a and b that we need to find the field value for
double Lerp(double a, double f_a, double b, double f_b, double x)
{
	double t = (x - a)/(b - a);
	return f_a + t*(f_b - f_a);
}

/*
 	This function sets our variables that we will use in main to represent the proper points 
	farPoint: the point that is not alligned by x value with any other point
	lowPoint: the lower of the two alligned points
	highPoint: the higher of the two alligned points
*/
void ConfigurePoints(const Triangle& triangle, Point &farPoint, Point &lowPoint, Point &highPoint, int x1, int x2)
{
	int x3 = 3-(x1+x2);				//3-(x1+x2) calculates the remaining point's index  
	highPoint.x = triangle.X[x1];
	highPoint.y = triangle.Y[x1];
	highPoint.z = triangle.Z[x1];
	highPoint.shading = triangle.shading[x1];

	lowPoint.x = triangle.X[x2];
	lowPoint.y = triangle.Y[x2];
	lowPoint.z = triangle.Z[x2];
	lowPoint.shading = triangle.shading[x2];

	farPoint.x = triangle.X[x3];
	farPoint.y = triangle.Y[x3];
	farPoint.z = triangle.Z[x3];
	farPoint.shading = triangle.shading[x3];

	for(int i=0; i<3; ++i)
	{
		highPoint.colors[i] = triangle.colors[x1][i];
		lowPoint.colors[i] = triangle.colors[x2][i];
		farPoint.colors[i] = triangle.colors[x3][i];
	}
}


/*
	This function figures out which two points share the same x value then of those two has the larger y value
	then uses that information to call ConfigurePoints
	x1: larger y value point
	x2: smaller y value point 
*/
void FindPoints(const Triangle &triangle, Point &farPoint, Point &lowPoint, Point &highPoint)
{	
	int x1, x2;              
	if(triangle.X[0] == triangle.X[1]) 
	{													
		if(triangle.Y[0] > triangle.Y[1])
		{
			x1=0;
			x2=1;
		}
		else
		{
			x1=1;
			x2=0;
		}
	}

	else if(triangle.X[0] == triangle.X[2]) 
	{
		if(triangle.Y[0] > triangle.Y[2])
		{
			x1=0;
			x2=2;
		}
		else
		{
			x1=2;
			x2=0;
		}
	}

	else
	{
		if(triangle.Y[2] > triangle.Y[1])
		{
			x1=2;
			x2=1;
		}
		else
		{
			x1=1;
			x2=2;
		}
	}

	ConfigurePoints(triangle, farPoint, lowPoint, highPoint, x1, x2);
	
}

/*
	This function performs the necessary checking for ranges that are out of bounds so we can
	calculate the amount of iterations of each inner loop we need accurately
	
	1C notes:
		- Modified to make modular by adding width and height parameters
*/
void CheckRangeBounds(double *xRange, double width, double height)
{
	
	if(xRange[0] < 0)													
		xRange[0] = 0;
	
	if(xRange[1] >= width)
		xRange[1] = width-1;
}

/*
	This function tells us the x and y ranges that we are interested in for each set of points that are passed to it

   xRange[0] == columnMin && xRange[1] == columnMax
	yRange[0] == rowMin && yRange[1] == rowMax
*/
void CalculateRanges(double *xRange, Point &highPoint, Point &lowPoint, Point &farPoint)
{

   if(lowPoint.x < farPoint.x)
   {
	   xRange[0] = ceil__mod(lowPoint.x);
	   xRange[1] = floor__mod(farPoint.x);
   }
   else
   {
	   xRange[0] = ceil__mod(farPoint.x);
	   xRange[1] = floor__mod(lowPoint.x);
   }
}

                                 /*PROJECT C SPECIFIC METHODS*/

/*
	This function returns true if the triangle has to points with the same x value
	likewise will return false if not

	If the triangle is not arbitrary there is less work to be done
*/

bool IsArbitrary(const Triangle& triangle)
{
	double x1, x2, x3;
	x1 = triangle.X[0];
	x2 = triangle.X[1];
	x3 = triangle.X[2];

	return !(x1==x2 || x1==x3 || x3==x2);

}

void FindEquationParameters(double leftX, double leftY, double rightX, double rightY, double *m, double *b)
{
	*m = (rightY - leftY)/(rightX - leftX);
	*b = leftY - ((*m)*leftX);

}

void Split(Triangle& triangle, Triangle* goingLeft, Triangle* goingRight)
{
	double x1, x2, x3, m, b, y;
	int leftIndex, rightIndex, middleIndex;
	x1 = triangle.X[0];
	x2 = triangle.X[1];
	x3 = triangle.X[2];

	//start finding the middle x
	if( (x2<x1 && x1<x3) || (x3<x1 && x1<x2) ) //if x1 is middle
	{
		if(x2<x3)
		{
			FindEquationParameters(x2, triangle.Y[1], x3, triangle.Y[2], &m, &b);
			leftIndex = 1;
			rightIndex = 2;
		}

		else
		{
			FindEquationParameters(x3, triangle.Y[2], x2, triangle.Y[1], &m, &b);
			leftIndex = 2;
			rightIndex = 1;
		}
		y = (m*x1)+b;
		middleIndex = 0;
	}

	else if( (x1<x2 && x2<x3) || (x3<x2 && x2<x1) ) //if x2 is middle
	{
		if(x1<x3)
		{
			FindEquationParameters(x1, triangle.Y[0], x3, triangle.Y[2], &m, &b);
			leftIndex = 0;
			rightIndex = 2;
		}
		else
		{
			FindEquationParameters(x3, triangle.Y[2], x1, triangle.Y[0], &m, &b);
			leftIndex = 2;
			rightIndex = 0;
		}
		y = (m*x2)+b;
		middleIndex = 1;
	}
	//if two x values are the same then this function wouldnt be called so this is just else
	else	//if x3 is middle
	{
		if(x1<x2)
		{
			FindEquationParameters(x1, triangle.Y[0], x2, triangle.Y[1], &m, &b);
			leftIndex = 0;
			rightIndex = 1;
		}
		else
		{
			FindEquationParameters(x2, triangle.Y[1], x1, triangle.Y[0], &m, &b);
			leftIndex = 1;
			rightIndex = 0;
		}
		y = (m*x3)+b;
		middleIndex = 2;
	}

//found middle x now split up triangles
	(*goingLeft).X[0] = triangle.X[leftIndex];
	(*goingLeft).Y[0] = triangle.Y[leftIndex];
	
	//get shading at this point
	(*goingLeft).shading[0] = triangle.shading[leftIndex];
	

	(*goingLeft).X[1] = triangle.X[middleIndex];
	(*goingLeft).Y[1] = triangle.Y[middleIndex];
	
	(*goingLeft).shading[1] = triangle.shading[middleIndex];

	(*goingLeft).X[2] = triangle.X[middleIndex];
	(*goingLeft).Y[2] = y;
	
	//lerp shading between points so no messy normal calculations later
	(*goingLeft).shading[2] = Lerp(triangle.X[leftIndex], triangle.shading[leftIndex], triangle.X[rightIndex], triangle.shading[rightIndex], triangle.X[middleIndex]);

	(*goingRight).X[0] = triangle.X[rightIndex];
	(*goingRight).Y[0] = triangle.Y[rightIndex];
	
	(*goingRight).shading[0] = triangle.shading[rightIndex];
	
	(*goingRight).X[1] = triangle.X[middleIndex];
	(*goingRight).Y[1] = triangle.Y[middleIndex];

	(*goingRight).shading[1] = triangle.shading[middleIndex];

	(*goingRight).X[2] = triangle.X[middleIndex];
	(*goingRight).Y[2] = y;

	(*goingRight).shading[2] = (*goingLeft).shading[2];
	
	//We made a new point so we need to lerp to figure out its color and z values
	double newPointColorLerp, newPointZLerp;
	
	//lerping new point's color field values
	for(int p=0; p<3; ++p) 			
	{	
		//lerp between left and right points for all triangles
		newPointColorLerp = Lerp(triangle.X[leftIndex], triangle.colors[leftIndex][p], triangle.X[rightIndex], triangle.colors[rightIndex][p], triangle.X[middleIndex]);
		
		(*goingRight).colors[2][p] = newPointColorLerp;	
		(*goingLeft).colors[2][p] = newPointColorLerp;
	}
	
	for(int p=0; p<3; ++p) //grab colors for goingLeft's left point and grab colors for goingRight's right point
	{
		(*goingLeft).colors[0][p] = triangle.colors[leftIndex][p];
		(*goingRight).colors[0][p] = triangle.colors[rightIndex][p];
	}

	for(int p=0; p<3; ++p) //similarly for above grab colors for both triangles middle point which they both have in common
	{
		(*goingLeft).colors[1][p] = triangle.colors[middleIndex][p];
		(*goingRight).colors[1][p] = triangle.colors[middleIndex][p];
	}
	
	//lerping new point's z field value
	newPointZLerp = Lerp(triangle.X[leftIndex], triangle.Z[leftIndex], triangle.X[rightIndex], triangle.Z[rightIndex], triangle.X[middleIndex]);
	
	(*goingRight).Z[0] = triangle.Z[rightIndex];
	(*goingRight).Z[1] = triangle.Z[middleIndex];
	(*goingRight).Z[2] = newPointZLerp;
	
	(*goingLeft).Z[0] = triangle.Z[leftIndex];
	(*goingLeft).Z[1] = triangle.Z[middleIndex];
	(*goingLeft).Z[2] = newPointZLerp;
	
}

void CalculateShading(Triangle &triangle, LightingParameters &lp, Camera &cam)
{
	for(int i=0; i<3; ++i)
	{
		double shadingValue = 0;
	//ambient
		
		shadingValue += lp.Ka;
		
 	//diffuse
		
		//dynamically allocated in normalize so delete after using
		double *vecN = Normalize(triangle.normals[i][0], triangle.normals[i][1], triangle.normals[i][2]);
		double vecL[3];
		
		vecL[0] = lp.lightDir[0];
		vecL[1] = lp.lightDir[1];
		vecL[2] = lp.lightDir[2];

		shadingValue += lp.Kd * max(0.0, DotProduct(vecL, vecN));	
	
	//specular
		
		double R[3], camVec[3], temp;
		temp = 2.0*DotProduct(vecL, vecN);
		
		for(int p=0; p<3; ++p)
		{
			vecN[p] *= temp;
		}
		
		for(int p=0; p<3; ++p)
		{
			R[p] = vecN[p] - vecL[p];
		}
	
		camVec[0] = cam.position[0] - triangle.X[i];
		camVec[1] = cam.position[1] - triangle.Y[i];
		camVec[2] = cam.position[2] - triangle.Z[i];

		Normalize(R);
		Normalize(camVec);
		
		shadingValue += lp.Ks * pow( max(0.0, DotProduct(R, camVec)), lp.alpha );
		
		delete vecN;
		triangle.shading[i] = shadingValue;
	}
	
}

int main()
{		
	
	std::vector<Triangle> triangles = GetTriangles(); 
   int npixels = HEIGHT*WIDTH*3;
   vtkImageData *image = NewImage(HEIGHT, WIDTH);
	unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);

	//put loop here for all four images for now just print one
	for(int c=0; c<1000; c++)
	{
		char imageName[128];
		sprintf(imageName,"frame%03d", c);
		Camera cam = GetCamera(c,1000);
		LightingParameters lp = GetLighting(cam);

      for (int i = 0 ; i < npixels ; i++)
          buffer[i] = 0;
	   
	   int zBuffSize = npixels/3;
	   double *zbuff = new double[zBuffSize];
	   for(int i=0; i<zBuffSize; ++i)
	   {
		   zbuff[i] = -1;
	   }
   
   
      Screen screen;
      screen.buffer = buffer;
      screen.width = WIDTH;
      screen.height = HEIGHT;
	   screen.zBuffer = zbuff;
		   
	   Matrix temp = Matrix::ComposeMatrices(cam.CameraTransform(), cam.ViewTransform());
	   Matrix m = Matrix::ComposeMatrices(temp, cam.DeviceTransform());
	   
	   
   
	   int size = triangles.size(); //upperBound-1 is the end of the buffer
      for(int p=0; p<size; ++p)                         //for every triangle
	   {
		   Triangle triangle = triangles[p];  
		   //instead of looping over triangles and applying transforms just do them for the individual triangles
		   //extract points from triangle while simultaneously converting them to 4d
		   CalculateShading(triangle, lp, cam);
			Point points[3];
		   Point tmpPnt;
		   for(int b=0; b<3; ++b)
		   {
			   points[b].x = triangle.X[b];
			   points[b].y = triangle.Y[b];
			   points[b].z = triangle.Z[b];
		   }
		   
		   //now transform point while you put the new point back in the copy of the triangle
		   for(int b=0; b<3; b++)
		   {
			   m.TransformPoint(&points[b], &tmpPnt);
			   triangle.X[b] = tmpPnt.x/tmpPnt.w;
			   triangle.Y[b] = tmpPnt.y/tmpPnt.w;
			   triangle.Z[b] = tmpPnt.z/tmpPnt.w;
		   }

	 // 	CalculateShading(triangle, lp, cam); 
		   int triSize;
		   Triangle** tris;
		   tris = new Triangle* [2];
		   tris[0] = nullptr;
		   tris[1] = nullptr;
   
		   Triangle goingLeft, goingRight;
		   if(IsArbitrary(triangle))
		   {
			   Split(triangle, &goingLeft, &goingRight);
			   triSize = 2;
			   tris[0] = &goingLeft;
			   tris[1] = &goingRight;
		   }
   
		   else
		   {
			   triSize = 1;
			   tris[0] = &triangle;
		   }
	   
		   for(int w=0; w<triSize; w++)
		   {	
			   if(tris[w] == nullptr)
				   exit(EXIT_FAILURE);
			   
			   Triangle& triangle2 = *(tris[w]);
			   
			   Point farPoint, lowPoint, highPoint;
			   double xRange[2];
				   
			   FindPoints(triangle2, farPoint, lowPoint, highPoint);
			   
			   CalculateRanges(xRange, highPoint, lowPoint, farPoint);
   
		      double m1, m2, b1, b2, y1, y2;            // y2: bottomEnd, y1: topEnd 
	         m1 = (farPoint.x > highPoint.x) ? 
		           (farPoint.y-highPoint.y)/(farPoint.x-highPoint.x) : (highPoint.y-farPoint.y)/(highPoint.x-farPoint.x);
	         
		      b1 = farPoint.y-(m1*farPoint.x); 
	         
		      m2 = (farPoint.x > lowPoint.x) ? 
			        (farPoint.y-lowPoint.y)/(farPoint.x-lowPoint.x) : (lowPoint.y-farPoint.y)/(lowPoint.x-farPoint.x);
	         
		      b2 = farPoint.y-(m2*farPoint.x); 
	         
		      CheckRangeBounds(xRange, WIDTH, HEIGHT); 		//check for bounds that are outside of the screen buffer
			   
		      for(int i=xRange[0]; i<=xRange[1]; ++i)              //for every scanline
	         {
	   	      y1 = floor__mod(b1+m1*i);
	   	      y2 = ceil__mod(b2+m2*i);
				   
				   /*
					   Making y1 and y2 Point objects makes operations easier
					   because we need to lerp their field values to find field
					   values for each pixel in each scanline for whom y1 and y2
					   are the top and bottom ends respectively
				   */
				   Point y1Point, y2Point;
				   y1Point.x = i;
				   y1Point.y = b1+m1*i;
				   y2Point.x = i;
				   y2Point.y = b2+m2*i;
				   
				   //find out if going left or going right because this will alter lerping operation
					//we are using the precalculated shading for each individual vertex to lerp the shading for the scanline end points
				   if(farPoint.x < lowPoint.x) //for going left 
				   {
					   for(int k=0; k<3; ++k)
					   {
						   y1Point.colors[k] = Lerp(farPoint.x, farPoint.colors[k], highPoint.x, highPoint.colors[k], i);
						   y2Point.colors[k] = Lerp(farPoint.x, farPoint.colors[k], lowPoint.x, lowPoint.colors[k], i);
					   }
					   y1Point.z = Lerp(farPoint.x, farPoint.z, highPoint.x, highPoint.z, y1Point.x);	
					   y1Point.shading = Lerp(farPoint.x, farPoint.shading, highPoint.x, highPoint.shading, y1Point.x);

					   y2Point.z = Lerp(farPoint.x, farPoint.z, lowPoint.x, lowPoint.z, y2Point.x);
					   y2Point.shading = Lerp(farPoint.x, farPoint.shading, lowPoint.x, lowPoint.shading, y2Point.x);
				   }
				   else //for going right
				   {
					   for(int k=0; k<3; ++k)
					   {
						   y1Point.colors[k] = Lerp(highPoint.x, highPoint.colors[k], farPoint.x, farPoint.colors[k], i);
						   y2Point.colors[k] = Lerp(lowPoint.x, lowPoint.colors[k], farPoint.x, farPoint.colors[k], i);
					   }
				   
					   y1Point.z = Lerp(highPoint.x, highPoint.z, farPoint.x, farPoint.z, y1Point.x);
					   y1Point.shading = Lerp(highPoint.x, highPoint.shading, farPoint.x, farPoint.shading, y1Point.x);

					   y2Point.z = Lerp(lowPoint.x, lowPoint.z, farPoint.x, farPoint.z, y2Point.x);
					   y2Point.shading = Lerp(lowPoint.x, lowPoint.shading, farPoint.x, farPoint.shading, y2Point.x);
				   }
	   	      
				   for(int j=(int)y2; j<=y1; ++j)            //for every pixel in the scanline
		         {
			         int index = 3*(j*WIDTH+i); 
					   double zLerp = Lerp(y2Point.y, y2Point.z, y1Point.y, y1Point.z, j);

					   //find out if pixel is even visable BEFORE doing color calculations
	   		      
					   if(index>=0 && index<npixels && screen.zBuffer[index/3] < zLerp)  //if the pixel is within the bound of the screen buffer then write colorss
				      {	
						   screen.zBuffer[index/3] = zLerp;
						   for(int k=0; k<3; ++k)
						   {	
								//we are lerping the shading value for each pixel based the shading values for the two end points of the scanline
								double shadeLerp = Lerp(y2Point.y, y2Point.shading, y1Point.y, y1Point.shading, j);
							   screen.buffer[index + k] = ceil__mod(255.0*min(1.0, shadeLerp*Lerp(y2Point.y, y2Point.colors[k], y1Point.y, y1Point.colors[k], j)));
						   }
						   
		   	      }		
			      }
		      }
		   }
	   }
      WriteImage(image, imageName);
	}
}
