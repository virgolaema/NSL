#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>

//chisquared, binned approach
double chisquared (int , double * , double * , double *);

//computing the variance in a MC method, using blocks
double variance_blocks (int, double *);

//distanza, prende in ingresso x1,x2,y1,y2
double distance(double, double, double, double);

//counts input values from file
int count_file (const char*);

//binary search for int
int binarySearch(int [], int, int, int);

//mergesort
void merge(int [],int,int,int);
void merge_sort(int[],int, int);

//resize a vector
void resize (double **, int, int );

//mean of a vector
double mean (int, double[]);

//variance of a vector
double variance (int, double[]);

//swap array values
void swapArray(double[], int , int );

//random numbers in interval a b
double randUnif (double, double);

//maxa value in a vector
double findMax (double[], int );

//boxmuller RIVEDERE
void boxMuller (double, double, double*);

//generate numbers from gaussian RIVEDERE
double* generarandgauss(double, double, int );

//ordina il vettore in modo decrescente
void ordinadecr (double[], int );