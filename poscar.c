#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<string.h>

//////////////////////////////////////////////////////////////////////////
//#define DEBUG

#define	LATTICE_CONSTANT	3.171512202

//the triple product of the basis vectors must be positive

#define A1	sqrt(0.5),	-sqrt(1.5),	0.0
#define A2	sqrt(0.5),	sqrt(1.5),	0.0
#define A3	0.0,		0.0,		sqrt(0.75)

#define NC	3 //number of atoms in origin cell
#define NX 	3
#define NY 	3
#define NZ 	10

#define NT	NX*NY*NZ//number of lattice points in super cell

#define	C1	0.0,		0.0,		0.0
#define	C2	2.0/3.0,	1.0/3.0,	1.0/3.0
#define	C3	1.0/3.0,	2.0/3.0,	2.0/3.0

#define VACUUM	10
#define FIX		0.49


///////////////////////////////////////////////////////////////////////////
double A123[3][3] = { { A1 },{ A2},{ A3} };
double C123[NC][3] = { { C1 },{ C2},{ C3} };

double cross(double a[], double b[])
{
	return	a[0] * b[0] + \
		a[1] * b[1] + \
		a[2] * b[2];
}
void vectors()
{
	double V[3];
    V[0]=A123[0][1]*A123[1][2]-A123[1][1]*A123[0][2];
    V[1]=A123[2][0]*A123[0][2]-A123[0][0]*A123[2][2];
    V[2]=A123[0][0]*A123[1][1]-A123[1][0]*A123[0][1];
    if(cross(V,A123[2])>0.0)printf("vectors right\n");
    else printf("the triple product of the basis vectors must be positive\n");
    return 0;
}



double length(double a[])
{
	return sqrt(cross(a, a));
}


void info()
{
	printf("the length (in angstrom) of the bulk super cell is\n\t%f\t%f\t%f\n",\
			length(A123[0])*NX*LATTICE_CONSTANT,\
			length(A123[1])*NY*LATTICE_CONSTANT,\
            length(A123[2])*NZ*LATTICE_CONSTANT);
    printf("the length (in angstrom) of the slab super cell is\n\t%f\t%f\t%f\n",\
			length(A123[0])*NX*LATTICE_CONSTANT,\
			length(A123[1])*NY*LATTICE_CONSTANT,\
            length(A123[2])*NZ*LATTICE_CONSTANT+VACUUM);
	return 0;
}

void poscar0()
{
	FILE *p = fopen("POSCAR0", "w");
	int s;

	fprintf(p, "POSCAR0\n");
	fprintf(p, "%20.15f\n",LATTICE_CONSTANT);
	fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[0][0],A123[0][1],A123[0][2]);
    fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[1][0],A123[1][1],A123[1][2]);
	fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[2][0],A123[2][1],A123[2][2]);
	fprintf(p, "%d\n", NC);
	fprintf(p, "Selective dynamics\n\tDirect\n");		
	for (s = 0; s < NC; ++s)
	{
		fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\tT\tT\tT\n", C123[s][0], C123[s][1], C123[s][2]);
	}

	fclose(p);
	return 0;						
}

void poscar_bulk(int x,int y,int z,char *fix)
{
	char filename[50]="POSCAR_BULK_\0";
    //puts(filename);
    strcat(filename,fix);

	
    
	FILE *p = fopen(filename, "w");
	int s,i,j,k;
   
    int f=0;

	fprintf(p, "POSCAR %dx%dx%d\n",x,y,z);
	fprintf(p, "%20.15f\n",LATTICE_CONSTANT);
	fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[0][0]*NX,A123[0][1]*NX,A123[0][2]*NX);
    fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[1][0]*NY,A123[1][1]*NY,A123[1][2]*NY);
	fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[2][0]*NZ,A123[2][1]*NZ,A123[2][2]*NZ);
	fprintf(p, "%d\n", NT*NC);
	fprintf(p, "Selective dynamics\n\tDirect\n");		
	for (i = 0; i < x; ++i)
	{
    for (j = 0; j < y; ++j)
	{
    for (k = 0; k < z; ++k)
	{
    for(s=0;s<NC;++s)
    {
		fprintf(p, "%20.15f\t%20.15f\t%20.15f\t%s\t%s\t%s\n",\
         (C123[s][0]+i)/NX, (C123[s][1]+j)/NY, (C123[s][2]+k)/NZ,\
         fix,fix,fix);
       
        ++f;
#ifdef DEBUG         
        printf("%d\t%d\t%d\t%d\t%d\n",i,j,k,s,f);
#endif
    }}}}

	if(f==NT*NC)printf("%d atoms created successfully\n",f);
	fclose(p);
	return 0;
}

void poscar_slab(int x,int y,int z,double vacuum,double fix)
{
	char filename[50]="POSCAR_SLAB\0";
	char T='T',F='F';
    //puts(filename);
    //strcat(filename,ftoa(vacuum));

	
    
	FILE *p = fopen(filename, "w");
	int s,i,j,k;
   
    int f=0,g=0;

	fprintf(p, "POSCAR %dx%dx%d %2.2f %2.2f\n",x,y,z,vacuum,fix);
	fprintf(p, "%20.15f\n",LATTICE_CONSTANT);
	fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[0][0]*NX,A123[0][1]*NX,A123[0][2]*NX);
    fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[1][0]*NY,A123[1][1]*NY,A123[1][2]*NY);
	fprintf(p, "\t%20.15f\t%20.15f\t%20.15f\n",A123[2][0]*NZ,A123[2][1]*NZ,A123[2][2]*NZ+vacuum/LATTICE_CONSTANT);
	fprintf(p, "%d\n", NT*NC);
	fprintf(p, "Selective dynamics\n\tDirect\n");		
	for (i = 0; i < x; ++i)
	{
    for (j = 0; j < y; ++j)
	{
    for (k = 0; k < z; ++k)
	{
    for(s=0;s<NC;++s)
    {
		if((C123[s][2]+k)/NZ>fix)
			fprintf(p, "%20.15f\t%20.15f\t%20.15f\t%c\t%c\t%c\n",\
					(C123[s][0]+i)/NX, (C123[s][1]+j)/NY, (C123[s][2]+k)/(NZ+vacuum/(LATTICE_CONSTANT*A123[2][2])),\
					T,T,T);
        else{
			fprintf(p, "%20.15f\t%20.15f\t%20.15f\t%c\t%c\t%c\n",\
					(C123[s][0]+i)/NX, (C123[s][1]+j)/NY, (C123[s][2]+k)/(NZ+vacuum/(LATTICE_CONSTANT*A123[2][2])),\
					F,F,F);
       
        ++f;}
        ++g;
#ifdef DEBUG         
        printf("%d\t%d\t%d\t%d\t%d\n",i,j,k,s,f);
#endif
    }}}}

	printf("%2.2f A vacuum added with %d buttom atoms (%2.2f%%) fixed\n",vacuum,f,(double)f/g*100);
	fclose(p);
	return 0;
}

int main()
{
	vectors();
	poscar0();
    poscar_bulk(NX,NY,NZ,"T");
    poscar_bulk(NX,NY,NZ,"F");
    poscar_slab(NX,NY,NZ,VACUUM,FIX);
	info();
    












	system("pause");
	return 0;
}
