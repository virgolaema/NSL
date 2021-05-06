#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include "random.h"

//factorial function
double factorial (double);

//mclaurin series of the cos, for small x   
double cos_ml(double);

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

//CLASSES 

//position

#ifndef __POSITION_H__
#define __POSITION_H__

class Position{
	public:
	Position();
	Position(double x, double y, double z);
	Position (Position&); //copy constructor
	//~Position();
	
	void copyPos (Position&);
    void setX(double x) {m_x = x;};
	void setY(double y) {m_y = y;};	
	void setZ(double z) {m_z = z;};
	double R() const;
	double Phi() const;
	double Theta() const;
	double Rho() const;
	double Distance(const Position&) const;
	double getX () const {return m_x;}
	double getY () const {return m_y;}
	double getZ () const {return m_z;}

	protected:
		double m_x,m_y,m_z;
};

#endif

