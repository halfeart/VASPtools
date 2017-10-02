#define _CRT_SECURE_NO_WARNINGS
//#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define LC 3.171512202 //simple box lattice constant seakmc:1 WH_Wang:3.14 WHe_Juslin:3.165

#define XX 2.0  //for seakmc always  > 0 no units
#define YY 2.99999999999999
#define ZZ 2.0

#define LL 8

#define DT 1.0e-10

#define NAMAX 10000000

#define A1	1,	2,	1
#define A2	-1,	0,	1
#define A3	1,	-1,	1

double A123[3][3] = { { A1 },{ A2},{ A3} };//xyz basis vectors//orthogonal

//[line][column] hang lie,each line for a vector,must be orthogonal

double LX, LY, LZ;


//#define BASIS1 "0 0 0"  
//#define BASIS2 "0.5 0.5 0.5"  //for BCC

double BASIS0[3] = { 0.0,0.0,0.0 };
double BASIS1[3] = { 0.5,0.5,0.5 };

double ATOM[NAMAX][3];
int NA = 0;

double cross(double a[], double b[])
{
	return	a[0] * b[0] + \
		a[1] * b[1] + \
		a[2] * b[2];
}

double length(double a[])
{
	return sqrt(cross(a, a));
}

int init()
{
	int i, j;
	double t = 0;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf("%f ", A123[i][j]);
		}printf("\n");
	}

	if (fabs(t = cross(A123[0], A123[1])) < DT&&\
		fabs(t = cross(A123[0], A123[2])) < DT&&\
		fabs(t = cross(A123[1], A123[2])) < DT)
		puts("vectors orthogonal.");
	else
		puts("wrong vectors!");

	LX = XX*length(A123[0]);
	LY = YY*length(A123[1]);
	LZ = ZZ*length(A123[2]);

	return 0;
}

int reduce()
{
	int i, j, k, s = 0;
	for (i = -LL * LX; i < LL * LX; i++)
		for (j = -LL * LY; j < LL * LY; j++)
			for (k = -LL * LZ; k < LL * LZ; k++)
			{
				ATOM[s][0] = i;
				ATOM[s][1] = j;
				ATOM[s][2] = k;
				s++;
				ATOM[s][0] = i + BASIS1[0];
				ATOM[s][1] = j + BASIS1[1];
				ATOM[s][2] = k + BASIS1[2];
				s++;
			}NA = s;
	/*
	for(s=0;s<NA;s++)
	printf("%f\t%f\t%f\n",ATOM[s][0],ATOM[s][1],ATOM[s][2]);
	*/
	return 0;
}


int lattice(int N)
{
	int s;
	double x, y, z;
	double l = LX*LY*LZ;
	int i;
	for (s = 0; s < NA; ++s)
	{
		x = ATOM[s][0] * A123[0][0] + ATOM[s][1] * A123[0][1] + ATOM[s][2] * A123[0][2];
		y = ATOM[s][0] * A123[1][0] + ATOM[s][1] * A123[1][1] + ATOM[s][2] * A123[1][2];
		z = ATOM[s][0] * A123[2][0] + ATOM[s][1] * A123[2][1] + ATOM[s][2] * A123[2][2];
		ATOM[s][0] = x / length(A123[0]);
		ATOM[s][1] = y / length(A123[1]);
		ATOM[s][2] = z / length(A123[2]);
	}
#ifdef DEBUG
	for (s = 0; s < NA; s++)
	{
		if (length(ATOM[s]) < l)
		{
			l = length(ATOM[s]);
			i = s;
		}
	}
	x = ATOM[i][0];
	y = ATOM[i][1];
	z = ATOM[i][2];
	puts("test");
#endif
	return 0;
}

int rotate(int N)
{

	int s = 0, i = 0;
	for (s = 0; s < N; s++)
	{
		if (ATOM[s][0] <= LX&&ATOM[s][1] <= LY&&ATOM[s][2] <= LZ)
		{
			if (ATOM[s][0] >= -DT&&ATOM[s][1] >= -DT&&ATOM[s][2] >= -DT)
			{
				ATOM[i][0] = ATOM[s][0];
				ATOM[i][1] = ATOM[s][1];
				ATOM[i][2] = ATOM[s][2];
				i++;
				//printf("%f\t%f\t%f\n",ATOM[s][0],ATOM[s][1],ATOM[s][2]);
			}
		}
	}NA = i; printf("%d\n", N);

	return 0;
}

int output(int N)
{
	int s;
	FILE *p = fopen("test.xyz", "w");
	fprintf(p, "%d\n\n", N);
	for (s = 0; s < N; ++s)
	{
		fprintf(p, "%f\t%f\t%f\t%d\n", ATOM[s][0], ATOM[s][1], ATOM[s][2], s + 1);
	}
	fclose(p);
	return 0;
}

int output_poscar()
{
	int s;
	FILE *p = fopen("POSCAR", "w");

	fprintf(p, "BULK W 111\n");
	fprintf(p, "%20.15f\n", LC);
	/*
	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", length(A123[0]), 0.0, 0.0);
	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", 0.0, length(A123[1]), 0.0);
	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", 0.0, 0.0, length(A123[2]));
	*/
	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", LX, 0.0, 0.0);
	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", 0.0, LY, 0.0);
	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", 0.0, 0.0, LZ);


	fprintf(p, "%d\n", NA);
	fprintf(p, "Selective dynamics\nDirect\n");
	for (s = 0; s < NA; ++s)
	{
		fprintf(p, "%20.15f\t%20.15f\t%20.15f\tF\tF\tF\n", ATOM[s][0] / LX, ATOM[s][1] / LY, ATOM[s][2] / LZ);
	}
	fclose(p);
	return 0;
}

int output_poscar_slab(double vacuum, double fix) //fix fraction of orgin cell from zero
{
	int s;
	FILE *p = fopen("POSCAR_SLAB", "w");

	fprintf(p, "SLAB W 111\n");
	fprintf(p, "%20.15f\n", LC);

	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", LX, 0.0, 0.0);
	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", 0.0, LY, 0.0);
	fprintf(p, "%20.15f\t%20.15f\t%20.15f\n", 0.0, 0.0, LZ + vacuum / LC);


	fprintf(p, "%d\n", NA);
	fprintf(p, "Selective dynamics\nDirect\n");
	for (s = 0; s < NA; ++s)
	{
		if (ATOM[s][2] / LZ > fix)
			fprintf(p, "%20.15f\t%20.15f\t%20.15f\tT\tT\tT\n", ATOM[s][0] / LX, ATOM[s][1] / LY, ATOM[s][2] / (LZ + vacuum / LC));
		else
			fprintf(p, "%20.15f\t%20.15f\t%20.15f\tF\tF\tF\n", ATOM[s][0] / LX, ATOM[s][1] / LY, ATOM[s][2] / (LZ + vacuum / LC));
	}
	fclose(p);

	return 0;
}

int minimize(char poscar[])
{//reload a poscar and use some potentials to minimze the structure stable and periodical 
	FILE *p = fopen("POSCAR", "r");
	char s[1000];
	fgets(s, 1000, p);
	puts(s);
	return 0;
}

int main()
{
	//printf("hello world!\n");
	init();
	//printf("A:%d",A123[1][2]);
	reduce();
	//printf("%d\n\n",NA);
	lattice(NA);
	rotate(NA);
	//////////////////////////////////////////
	output(NA);
	output_poscar();
	output_poscar_slab(10.0, 1.0/3.0);

	minimize("POSCAR");

	return 0;
}
