#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#ifdef _WIN32
#include <conio.h>

#endif
double **TensorOfInertia;
double ** TensorDiag;
double *eigenvalues;
double **eigenvectors;

int* Types;

double TensorParameters[3];
double CONTACT_CUT = 2.5;

double* rgs; //starts with 2
double* rs;
double* conts;
double* rgwithout2s;



int N;
int Nrot;

struct Vector;





//old matrix eigenvals/eigenvecs code



#define ABS(x) (( (x)>0 )? (x) : -(x) )
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
#define MACHINE_ZERO 1.0e-35

/*	Gegeben a, berechnet Jacobi seine EW-e und die EV-en.
	Am Ende: d enthaelt die Eigenwerte und
				v die normierten Eigenvectoren als Spalten
	Indexbereich: 1..n */

/*Djacobi macht dasselbe wie Jacobi fuer double-Variable*/
void jacobi_double(double  **a, int n, double  d[], double  **v, int *nrot) {

	int j,iq,ip,i;
	double  tresh,theta,tau,t,sm,s,h,g,c,*b=NULL,*z=NULL;

//	void nrerror(char *);

/*
printf("jacobi = %lg %lg %lg \n",a[1][1],a[1][2],a[1][3]);
printf("jacobi = %lg %lg %lg \n",a[2][1],a[2][2],a[2][3]);
printf("jacobi = %lg %lg %lg \n",a[3][1],a[3][2],a[3][3]);
exit(1);
*/

	b=(double *)calloc(sizeof(double), n+1);
	z=(double *)calloc(sizeof(double), n+1);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1; ip<=n; ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1; i<=50; i++) { /* try  50 sweeps at most */
		sm=0.0;
		for (ip=1; ip<=n-1; ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += ABS(a[ip][iq]);
		}
		/*if(sm == 0.0) {  */
		/*The comparation of sm with zero wouldn't work, if the underflow is not */
		/*automatically set to zero. So we defined a number called MACHINE_ZERO: */
		if (sm <= MACHINE_ZERO ) {
			free(z);
			free(b);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*ABS(a[ip][iq]);
				if (i > 4 && (double )(ABS(d[ip])+g) == (double )ABS(d[ip])
					&& (double )(ABS(d[iq])+g) == (double )ABS(d[iq]))
					a[ip][iq]=0.0;
				else if (ABS(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double )(ABS(h)+g) == (double )ABS(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(ABS(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	printf("Too many iterations in routine jacobi");
}
#undef ROTATE

/* Sortiert die Eigenwerte im Array d[1..n] nach fallender Reihenfolge
   und Ordnet die  Entsprechenden Eigenvektoren in der Matrix v um */
void eigsort_in_decending_order(double d[], double **v, int n)
{
	int k,j,i;
	double p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] > p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}



//end old code





struct rgb{
    double r;       // percent
    double g;       // percent
    double b;       // percent
} ;
 struct hsv{
    double h;       // angle in degrees
    double s;       // percent
    double v;       // percent
} ;

hsv   rgb2hsv(rgb in);
rgb   hsv2rgb(hsv in);

hsv rgb2hsv(rgb in)
{
    hsv         out;
    double      min, max, delta;

    min = in.r < in.g ? in.r : in.g;
    min = min  < in.b ? min  : in.b;

    max = in.r > in.g ? in.r : in.g;
    max = max  > in.b ? max  : in.b;

    out.v = max;                                // v
    delta = max - min;
    if (delta < 0.00001)
    {
        out.s = 0;
        out.h = 0; // undefined, maybe nan?
        return out;
    }
    if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
        out.s = (delta / max);                  // s
    } else {
        // if max is 0, then r = g = b = 0              
            // s = 0, v is undefined
        out.s = 0.0;
        out.h = 0;                            // its now undefined
        return out;
    }
    if( in.r >= max )                           // > is bogus, just keeps compilor happy
        out.h = ( in.g - in.b ) / delta;        // between yellow & magenta
    else
    if( in.g >= max )
        out.h = 2.0 + ( in.b - in.r ) / delta;  // between cyan & yellow
    else
        out.h = 4.0 + ( in.r - in.g ) / delta;  // between magenta & cyan

    out.h *= 60.0;                              // degrees

    if( out.h < 0.0 )
        out.h += 360.0;

    return out;
}


rgb hsv2rgb(hsv in)
{
    double      hh, p, q, t, ff;
    long        i;
    rgb         out;

    if(in.s <= 0.0) {       // < is bogus, just shuts up warnings
        out.r = in.v;
        out.g = in.v;
        out.b = in.v;
        return out;
    }
    hh = in.h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.v * (1.0 - in.s);
    q = in.v * (1.0 - (in.s * ff));
    t = in.v * (1.0 - (in.s * (1.0 - ff)));

    switch(i) {
    case 0:
        out.r = in.v;
        out.g = t;
        out.b = p;
        break;
    case 1:
        out.r = q;
        out.g = in.v;
        out.b = p;
        break;
    case 2:
        out.r = p;
        out.g = in.v;
        out.b = t;
        break;

    case 3:
        out.r = p;
        out.g = q;
        out.b = in.v;
        break;
    case 4:
        out.r = t;
        out.g = p;
        out.b = in.v;
        break;
    case 5:
    default:
        out.r = in.v;
        out.g = p;
        out.b = q;
        break;
    }
    return out;     
}






struct Vector
{
public: double x[3];
       double L;
       Vector(double x1,double y1,double z1)
       {
                  x[0] = x1;
                  x[1] = y1;
                  x[2] = z1;

       //           L = sqrt(x1*x1 + y1*y1 + z1*z1);

       }
       Vector()
       {
       }


       Vector& operator = (const Vector &n)
       {
              x[0] = n.x[0];
              x[1] = n.x[1];
              x[2] = n.x[2];

 //                 L = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
       return *this;
       }
       bool operator == (const Vector &n)
       {
       if(x[0] == n.x[0] && x[1] == n.x[1] && x[2] == n.x[2])
       {
       return true;
       }
       else
       {
       return false;
       }
       }
       double len()
       {
       return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
       }
       double sqlen()
       {
           return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
       }

      /*  Vector getReal()
        {
            Vector tmp;
            for(int i = 0; i < 3; i++)
            {

            tmp.x[i]=x[i]-SIZE[i]*floor(x[i]/SIZE[i]);


            }

        return tmp;
        }
*/



};

Vector operator+(const Vector &a,const Vector &b)
       {

          Vector c(a.x[0]+b.x[0],a.x[1]+b.x[1],a.x[2]+b.x[2]);
          return c;
       }
Vector operator-(const Vector &a,const Vector &b)
       {
          Vector c(a.x[0]-b.x[0],a.x[1]-b.x[1],a.x[2]-b.x[2]);
          return c;
       }


/*
Vector getReal(const Vector &un)
        {
        Vector tmp;
        for(int i = 0; i < 3; i++)
        {

        tmp.x[i]=un.x[i]-SIZE[i]*floor(un.x[i]/SIZE[i]);


        }

        return tmp;
        }
*/
double operator*(const Vector &a,const Vector &b)
{
     double c = a.x[0]*b.x[0] +a.x[1]*b.x[1]+a.x[2]*b.x[2];
          return c;

}

Vector operator*(const double a,const Vector &b)
{
 Vector tmp;
 tmp = b;
 for(int i = 0; i < 3; i++)
        {

        tmp.x[i]=a*tmp.x[i];


        }
return tmp;
}


Vector* crd;




Vector coor[1000000];
Vector realCoor[1000000];
Vector endtoend[500000];
double boxSize[3] = {40,40,40};
int lcLength = 7;
int nLC;
int nAtoms;

	
	

double nearestImageR(const Vector &a,const Vector &b)
{
    Vector coor;
    for(int i=0; i < 3; i++)
    {
        coor.x[i] = a.x[i]-b.x[i];
        coor.x[i] = coor.x[i] - boxSize[i]*round(coor.x[i]/boxSize[i]); //PBC
    }
    return coor.len();
}
	
struct outp
{
double nemParam;
Vector direct;
};

outp localNematicParameter(Vector realCoor[1000000], double cutoff, Vector placement);

	
int main(int argc, char *argv[])
{
	
	
	char filename[100];
	FILE *fparam = fopen("params.txt","r");
	fscanf(fparam,"%s",&filename);
	fscanf(fparam,"%lf %lf %lf",&boxSize[0],&boxSize[1],&boxSize[2]);
	fscanf(fparam,"%i",&lcLength);
	FILE *fp = fopen(filename,"r");
	
	 nAtoms = 0;
	
	
	
	while(fscanf(fp,"%lf   %lf   %lf", &coor[nAtoms].x[0], &coor[nAtoms].x[1], &coor[nAtoms].x[2] )!=EOF)
	{
	//	coor[nAtoms].x[0] = (double)rand()/(double)boxSize[0];
	//	coor[nAtoms].x[1] = (double)rand()/(double)boxSize[0];
	//	coor[nAtoms].x[2] = (double)rand()/(double)boxSize[0];

		nAtoms++;
	}
	printf("read %i atoms\n",nAtoms);
	//PBC unwrapping
	for(int i = 0; i < nAtoms/lcLength; i++ )
	{
		realCoor[i*lcLength] = coor[i*lcLength];
		for(int j = i*lcLength+1; j < (i+1)*lcLength; j ++)
		{
			
			for(int d = 0; d < 3; d++)
			{
				
				if( fabs(coor[j].x[d] - realCoor[j-1].x[d]) > 5) 
				{
					//printf("correcting pbc....\n");
					if(coor[j].x[d] - realCoor[j-1].x[d] > 0)
					{
						realCoor[j].x[d] = coor[j].x[d] - boxSize[d];
					}
					else
					{
						realCoor[j].x[d] = coor[j].x[d] + boxSize[d];
						
					}
					if( fabs(realCoor[j].x[d] - realCoor[j-1].x[d]) > 5) 
					{
						printf("error unwrapping PBC\n");
					}
				}
				else
				{
					realCoor[j].x[d] = coor[j].x[d];
				}
			}				
			
		}
		endtoend[i] = (realCoor[i*lcLength] - realCoor[(i+1)*lcLength - 1]);
		if(endtoend[i].len() == 0) printf("error\n");
		endtoend[i] = (1.0/endtoend[i].len())*endtoend[i];
		for(int k = 0; k < 3; k++)
		{
			//if(endtoend[i].x[k] < 0) endtoend[i].x[k] = -endtoend[i].x[k];
		}
	}
	//end unwrapping
	nLC = nAtoms/lcLength;
//	Vector director = new Vector(0,0,0);
	printf("number of liquid crystal molecules %i\n",nLC);
	double **NematicTensor;
	NematicTensor = new double*[4];
	double *eigenvalues = new double[4];
	double **eigenvectors = new double*[4];
	for(int ss = 0; ss < 4; ss++)
	{
		eigenvectors[ss] = new double[4];
		NematicTensor[ss] = new double[4];
	}
	
	
	for(int xx = 0; xx < 4; xx++)
		{
			for(int yy = 0; yy < 4; yy++)
			{
				NematicTensor[xx][yy] = 0;
				
			}
			
		}
	for(int k = 0; k < nLC; k++)
	{
		for(int xx = 0; xx < 3; xx++)
		{
			for(int yy = 0; yy < 3; yy++)
			{
				if( xx == yy)
				{
					NematicTensor[xx+1][yy+1] -= 1.0/(3.0*(double)nLC);
				}
				NematicTensor[xx+1][yy+1] += endtoend[k].x[xx]*endtoend[k].x[yy]/(double)nLC;
			}
			
		}
		
	}
	
		for(int xx = 0; xx < 4; xx++)
		{
			for(int yy = 0; yy < 4; yy++)
			{
				printf("%lf ",NematicTensor[xx][yy]);
				
			}
			printf("\n");
		}
	
	
	int Nrot;
	 jacobi_double(NematicTensor, 3, eigenvalues,eigenvectors, &Nrot) ;
		eigsort_in_decending_order(eigenvalues,eigenvectors,3);
		printf("calculated tensor\n");
	
	printf("eigenvalues %lf %lf %lf\n" ,eigenvalues[1], eigenvalues[2],eigenvalues[3]);
	for(int m = 1; m < 4; m++)
	{
		printf("eigenvector %lf %lf %lf\n", eigenvectors[m][1],eigenvectors[m][2],eigenvectors[m][3]);
	}
	/*
	for(int k = 0; k < nLC; k +=200)
	{
		printf("sample endtoend %lf %lf %lf\n", endtoend[k].x[0], endtoend[k].x[1],endtoend[k].x[2]);
	}
	*/
	//largest eigenvalue eigenvector
	Vector director(eigenvectors[1][1],eigenvectors[2][1],eigenvectors[3][1]);
	
	
	
	//	printf("sample endtoend %lf %lf %lf\n", endtoend[0].x[0], endtoend[0].x[1],endtoend[0].x[2]);
//	printf("sample endtoend %lf %lf %lf\n", endtoend[4].x[0], endtoend[4].x[1],endtoend[4].x[2]);
//	printf("sample endtoend %lf %lf %lf\n", endtoend[7].x[0], endtoend[7].x[1],endtoend[7].x[2]);
//	printf("sample endtoend %lf %lf %lf\n", endtoend[20].x[0], endtoend[20].x[1],endtoend[20].x[2]);
	
	
	printf("%lf %lf %lf director\n",director.x[0],director.x[1],director.x[2]);
	
	double nemParam = 0; 
	for(int i = 0; i < nLC; i++)
	{
		nemParam += (3*(director*endtoend[i])*(director*endtoend[i])-1)/2;
	}
	printf("%lf all-cell nematic parameter\n",nemParam/(double)nLC);
	
	
	//local nematic parameters
	double cutoff = 9;
	double meanLocalParameter = 0;
	int cnt = 0;
	double minLocalParameter = 1;
	double mini,minj,mink;
	Vector field[50][50][50];
	int ci = 0,cj = 0,ck = 0;
	FILE *fhist = fopen("loc_hist.txt","w");
	for(double i = 0.1; i < boxSize[0]; i+=2*cutoff)
	{
		
		cj = 0;	
		for(double j = 0.1; j < boxSize[1]; j+=2*cutoff)
		{
			ck = 0;
			for(double k = 0.1; k < boxSize[2]; k+=2*cutoff)
			{
				
				Vector point(i,j,k);
				outp locNem = localNematicParameter(realCoor,cutoff,point);
				meanLocalParameter += locNem.nemParam;
				fprintf(fhist,"%lf\n",locNem.nemParam);
				if(locNem.nemParam < minLocalParameter)
				{
					mini = i; minj = j; mink = k;
					minLocalParameter = locNem.nemParam;				
				}
				cnt++;
				field[ci][cj][ck] = locNem.direct;
			ck++;
			}
		cj++;
		}
	ci++;
	}

	meanLocalParameter = meanLocalParameter/(double)cnt;
	
	printf("mean local nematic parameter with cutoff %lf is %lf\n",cutoff,meanLocalParameter);
	printf("minLocalParameter is %lf at %lf %lf %lf", minLocalParameter, mini,minj,mink);
	
	for(int k = 0; k < ck ; k++)
	{
	char str[80];
	sprintf(str,"sliceZ%i.txt",k);
	FILE* fk = fopen(str,"w");
	
	for(int i = 0; i < ci; i++)
	{
		for(int j = 0; j < cj; j++)
		{
			fprintf(fk,"%i %i %lf %lf \n",i,j, field[i][j][k].x[0], field[i][j][k].x[1] );
			
		}
		//printf("\n");	
	}
	fclose(fk);
	}
	
	for(int k = 0; k < ck ; k++)
	{
	char str[80];
	sprintf(str,"sliceX%i.txt",k);
	FILE* fk = fopen(str,"w");
	
	for(int i = 0; i < ci; i++)
	{
		for(int j = 0; j < cj; j++)
		{
			fprintf(fk,"%i %i %lf %lf \n",i,j, field[k][i][j].x[1], field[k][i][j].x[2] );
			
		}
		//printf("\n");	
	}
	fclose(fk);
	}
	
	for(int k = 0; k < ck ; k++)
	{
	char str[80];
	sprintf(str,"sliceY%i.txt",k);
	FILE* fk = fopen(str,"w");
	
	for(int i = 0; i < ci; i++)
	{
		for(int j = 0; j < cj; j++)
		{
			fprintf(fk,"%i %i %lf %lf \n",i,j, field[i][k][j].x[0], field[i][k][j].x[2] );
			
		}
		//printf("\n");	
	}
	fclose(fk);
	}
	
	
	
	
	
	
	return 0;
}


outp localNematicParameter(Vector realCoor[1000000], double cutoff, Vector placement)
{
	int nLC = nAtoms/lcLength;
	bool LCNeighbors[50000];
	int NLCNeighbors = 0;
	Vector endtoend[50000];
	int counterLC = 0;
	for(int i = 0; i < nLC; i++ )
	{
		LCNeighbors[i] = false;
		for(int j = i*lcLength; j < (i+1)*lcLength; j++)
		{
		if( nearestImageR(realCoor[j], placement) < cutoff)
		{
			LCNeighbors[i] = true;
			counterLC++;
		}
		}
		endtoend[i] = (realCoor[i*lcLength] - realCoor[(i+1)*lcLength - 1]);
		if(endtoend[i].len() == 0) printf("error\n");
		endtoend[i] = (1.0/endtoend[i].len())*endtoend[i];
	}
	
	if(counterLC < 5)
	{
	outp def;
	def.direct.x[0] = 0;
	def.direct.x[1] = 0;
	def.direct.x[2] = 0;
	def.nemParam = 0;
	return def;
	}
	
	
	double **NematicTensor;
	NematicTensor = new double*[4];
	double *eigenvalues = new double[4];
	double **eigenvectors = new double*[4];
	for(int ss = 0; ss < 4; ss++)
	{
		eigenvectors[ss] = new double[4];
		NematicTensor[ss] = new double[4];
	}
	
	
	for(int xx = 0; xx < 4; xx++)
		{
			for(int yy = 0; yy < 4; yy++)
			{
				NematicTensor[xx][yy] = 0;
				
			}
			
		}
	for(int k = 0; k < nLC; k++)
	{
		if(LCNeighbors[k] == true)
		{
		NLCNeighbors++;
		for(int xx = 0; xx < 3; xx++)
		{
			for(int yy = 0; yy < 3; yy++)
			{
				if( xx == yy)
				{
					NematicTensor[xx+1][yy+1] -= 1.0/3.0;
				}
				NematicTensor[xx+1][yy+1] += endtoend[k].x[xx]*endtoend[k].x[yy];
			}
			
		}
		}
	}
	
		for(int xx = 0; xx < 4; xx++)
		{
			for(int yy = 0; yy < 4; yy++)
			{
				NematicTensor[xx][yy] = NematicTensor[xx][yy]/NLCNeighbors;
		//		printf("%lf ",NematicTensor[xx][yy]);
				
			}
		//	printf("\n");
		}
	
	
	int Nrot;
	 jacobi_double(NematicTensor, 3, eigenvalues,eigenvectors, &Nrot) ;
		eigsort_in_decending_order(eigenvalues,eigenvectors,3);
	//	printf("calculated tensor\n");
	
//	printf("eigenvalues %lf %lf %lf\n" ,eigenvalues[1], eigenvalues[2],eigenvalues[3]);
	for(int m = 1; m < 4; m++)
	{
	//	printf("eigenvector %lf %lf %lf\n", eigenvectors[m][1],eigenvectors[m][2],eigenvectors[m][3]);
	}
	/*
	for(int k = 0; k < nLC; k +=200)
	{
		printf("sample endtoend %lf %lf %lf\n", endtoend[k].x[0], endtoend[k].x[1],endtoend[k].x[2]);
	}
	*/
	//largest eigenvalue eigenvector
	Vector director(eigenvectors[1][1],eigenvectors[2][1],eigenvectors[3][1]);
	
	
	
	//	printf("sample endtoend %lf %lf %lf\n", endtoend[0].x[0], endtoend[0].x[1],endtoend[0].x[2]);
//	printf("sample endtoend %lf %lf %lf\n", endtoend[4].x[0], endtoend[4].x[1],endtoend[4].x[2]);
//	printf("sample endtoend %lf %lf %lf\n", endtoend[7].x[0], endtoend[7].x[1],endtoend[7].x[2]);
//	printf("sample endtoend %lf %lf %lf\n", endtoend[20].x[0], endtoend[20].x[1],endtoend[20].x[2]);
	
	
	//printf("%lf %lf %lf director\n",director.x[0],director.x[1],director.x[2]);
	
	double nemParam = 0; 
	for(int i = 0; i < nLC; i++)
	{
		if(LCNeighbors[i] == true)
		{
		nemParam += (3*(director*endtoend[i])*(director*endtoend[i])-1)/2;
		}
	}
//	printf("%lf all-cell nematic parameter\n",nemParam/(double)nLC);
	

	for(int ss = 0; ss < 4; ss++)
	{
		delete[] eigenvectors[ss] ;
		delete[] NematicTensor[ss] ;
	}
	delete[] NematicTensor;
	
	delete[] eigenvalues ;
	delete[] eigenvectors ;
	outp xx;
	xx.direct = director;
	if(director.len() < 0.5) printf("error");
	xx.nemParam = nemParam/(double)NLCNeighbors;
	return xx;
	
	
}


