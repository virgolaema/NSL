#include "lib.h"

using namespace std;

//factorial function
double factorial (double x){
    int fact = 1;
    int appo = x;
    if (x==0){appo=1;}
    while (appo!=1){
        fact *= appo;
        appo--;
    }
    return fact;
}

//mclaurin series of the cos, for small x   
double cos_ml(double x){
    double sum = 0;
    for (int i=0; i<3; i++){ //first terms of the series
        sum +=  pow(-1,i)/factorial(2*i)*pow(x,2*i);
    }
    return sum;
}

//chisquared, binned approach
double chisquared (int N, double * nu, double * values, double * variances){
    double chisq = 0;
    for(int i = 0; i < N; i++){
        chisq += pow(nu[i] - values[i],2)/variances[i];
    }
    return chisq;
}

//computing the variance in a MC method, using blocks
double variance_blocks (int N_blocks, double * A){
    double A2_mean = 0;
    for (int i=0; i<N_blocks; i++)  A2_mean += pow(A[i],2);
    A2_mean /= N_blocks;
    double A_mean = mean (N_blocks, A);
    return 1./(N_blocks-1)*(A2_mean-A_mean*A_mean);
}

//distance between 2 points in 2d space
double distance (double x1, double x2, double y1, double y2){
    double d, xx, yy;
    xx  =  pow(x1-x2,2);
    yy  =  pow(y1-y2,2);
    d  =  pow(xx + yy, 0.5);
    return d;
}

//counts input values from file
int count_file (const char* filename){
	ifstream in;
	int c = 0;
	double appo;
	in.open(filename);
    	if(in.fail()){
        cout  <<  endl  << "Can't open file!!!!"  <<  endl;
        return 1;
	} 

    in>>appo;
    while(!in.eof()){
	   c++;
	   in>>appo;
    }

    in.close();
    return c;
}

//binary search for int
int binarySearch(int v[], int value, int low, int high) {
     int mid;
     if (high < low){
	cout << "Impossible search!!!" << endl;
         return -1; 
	}

     mid  =  low + ((high - low) / 2); 
     if (v[mid] > value)
         return binarySearch(v, value, low, mid-1);
     else if (v[mid] < value)
         return binarySearch(v, value, mid+1, high);
     else
         return mid; 
 }


//mergesort
void merge(int a[],int low,int mid,int high){

    int h,i,j,k;
    int b[50]; 
    h = low; 
    i = low; 
    j = mid+1; 

    while((h<= mid)&&(j<= high)){ 
        if(a[h]<= a[j]){
            b[i] = a[h];
            h++; 
        }
        else{
            b[i] = a[j];
            j++;
        }
        i++;
    }

    if(h>mid){ 
        for(k = j;k<= high;k++){
            b[i] = a[k];
            i++;
        }
    }
    else{
        for(k = h;k<= mid;k++){ 
            b[i] = a[k];
            i++;
        }
    }
    for(k = low;k<= high;k++) a[k] = b[k]; 
}

void merge_sort(int a[],int low,int high){
    int mid;
    if(low<high) {
        mid  =  low + (high-low)/2; 
        merge_sort(a,low,mid);
        merge_sort(a,mid+1,high);
        merge(a,low,mid,high);
    }
}

//resize a vector
void resize (double **v, int oldDim, int newDim){
	double*pappo;

	if(newDim<oldDim){
		cout << "Set a larger dimension than previous one!" << endl;
		return;
	}

	pappo =  new double [newDim];

	for (int i = 0; i<oldDim; i++){
		pappo[i] = (*v)[i];
	}

	delete [] *v;
	*v = pappo;
    delete [] pappo;
return;
}

//mean of a vector
double mean(int a, double vett[]){
    double c,tot = 0;

    for (int i = 0; i<a; i++){
        tot =  tot+ vett[i];
    }

    c = tot/a;
    return c;
}


//variance of a vector
double variance (int a, double vett[]){
    double tot2 = 0;
    double c,d;

    for (int i = 0; i<a; i++){
        tot2 =  tot2+ vett[i]*vett[i];
    }
    d = mean(a, vett);
    c =  tot2/a - d*d;
    return c;
}

//swap array values
void swapArray(double v[], int pos1, int pos2){
    double appo;
    appo  =  v[pos1];
    v[pos1]  =  v[pos2];
    v[pos2]  =  appo;
}

//random numbers in interval a b
double randUnif (double a ,double b){
	double appo;
	appo = rand();
	appo = appo/RAND_MAX;
	appo = appo*(a-b)+b;
    return appo;
}


//max value in a vector
double findMax (double vett[], int dim) {
    int i;
    double max = 0;
    for (i = 0; i<dim; i++){
        if (vett[i]>max){
            max = vett[i]; 
        }
    }
    return max;
}

//boxmuller 
double boxmuller (double mean, double stdDev){
    Random rnd;
    rnd.Initialize(rnd);
    double appo1 =  rnd.Rannyu();
    double appo2 =  rnd.Rannyu();
    rnd.SaveSeed();
    return  stdDev*(sqrt(-2*log(appo1))*sin(2*(3.14)*appo2))+mean;
}

//ordina il vettore in ordine decrescente
void ordinadecr (double vett[], int dim){
  
    int min = 0;
    double temp;
    
    for(int i = 0; i<dim-1; i++){
        
        min  =  i;
        
        for(int j = i+1; j<dim; j++)
            if(vett[j] > vett[min])
                min  =  j;
        
        temp = vett[min];
        vett[min] = vett[i];
        vett[i] = temp;
    }
    
}

//CLASSES

//Position

Position::Position(){
    m_x=0;
    m_y=0;
    m_z=0;
}

Position::Position(double x, double y, double z){
    m_x=x;
    m_y=y;
    m_z=z;
}

Position::Position (Position &p){
    m_x = p.getX();
    m_y = p.getY();
    m_z = p.getZ();
}    

Position::~Position(){}

void Position::copyPos (Position &p){
    m_x = p.getX();
    m_y = p.getY();
    m_z = p.getZ();
}

void Position :: setCoor (double x, double y, double z){
    m_x = x;
    m_z = z;
    m_y = y;
}

double Position::Distance (const Position &p) const{return sqrt(pow(m_x-p.getX(),2) + pow(m_y-p.getY(),2) + pow(m_z-p.getZ(),2));}
double Position::R() const {return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);}
double Position::Phi () const {	return atan2(m_y,m_x);}
double Position::Theta () const {return acos(m_z/R());}
double Position::Rho () const {return sqrt(m_x*m_x+m_y*m_y);}

