#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
using namespace std;

void get_peaks(double x[1329], double a[1329]){

    double y_prev = a[0];

    for(int i = 1; i < 1329; i++){
        if(y_prev < 0 and a[i] > 0)
            cout << "start between " << x[i-1] << " and " << x[i] << endl;
        if(y_prev > 0 and a[i] < 0)
            cout << "end between " << x[i-1] << " and " << x[i] << endl;
        y_prev = a[i];
    }

}

/*
double bisection(double a, double b){

    cout << "Bisection Algorithm" << endl;

    int n = 1;

    double fa = f_x(a);

    cout << "Iteration\t" << "Value of x\t" << "Value of f(x)\t" << "Absolute Error\t" << "Relative Error" <<  endl;

    while( n <= 100 ){
        double p = a + (b-a)/2.0;

        cout << n << "\t";

        printf("%15f\t", p);

        double fp = f_x(p);

        double abs_err = abs(truth - p);

        printf("%15f\t%15f\t%15f\n", fp, abs_err, abs_err/abs(truth));

        if(fp == 0 || (b-a)/2.0 < TOL){
            return p;
        }

        if (fa*fp > 0){
            a = p;
            fa = fp;
        }

        else{
            b = p;
        }

        n += 1;
     }

     return NULL;
}
*/

void sg_filter(int nFilterPoints, int nPasses, double x[1329], double a[1329]){

    double c5[5] = {-3, 12, 17, 12, -3};
    double c11[11] = {-36, 9, 44, 69, 84, 89, 84, 69, 44, 9, -36};
    double c17[17] = {-21, -6, 7, 18, 27, 34, 39, 42, 43, 42, 39, 34, 27, 18, 7, -6, -21};

    for(int p = 1; p <= nPasses; p++)
        for(int i = (nFilterPoints-1)/2; i < 1329-(nFilterPoints-1)/2; i++){
            double sumX = 0.0;
            double sumA = 0.0;
            for(int j = 0; j < nFilterPoints; j++){
                
                if(nFilterPoints == 5){
                    sumA += c5[j]*a[i+1];
                }
            
                else if(nFilterPoints == 11){
                    sumA += c11[j]*a[i+1];
                }

                else if(nFilterPoints == 17){
                    sumA += c17[j]*a[i+1];
                }
                //sumA += a[j];
            }
            //if (p == nPasses)
                //cout << x[i] << "\t" << sumA/nFilterPoints << endl;
                a[i] = sumA/nFilterPoints;    
    }

}


void boxcar_filter(int nFilterPoints, int nPasses, double x[1329], double a[1329]){
    
   for(int m = 0; m < nPasses; m++)
    for(int i = (nFilterPoints-1)/2; i < 1329-(nFilterPoints-1)/2; i++){
        double sumX = 0.0;
        double sumA = 0.0;
        for(int j = i - (nFilterPoints-1)/2; j <= i + (nFilterPoints-1)/2; j++){
            sumA += a[j];
        }
        //if (m == 4)
            //cout << x[i] << "\t" << sumA/nFilterPoints << endl;
        a[i] = sumA/nFilterPoints;

    }
    

}

void getCoefficients(int n, double x[1329], double a[1329], double b[1329], double c[1329], double d[1329]){

    double h[n+1];
    double alpha[n+1];
    double l[n+1];
    double mu[n+1];
    double z[n+1];
    
    for(int i = 0; i < n; i++){
        h[i] = x[i+1] - x[i];
    }

    for(int i = 1; i < n; i++){
        alpha[i] = (3.0/h[i])*(a[i+1]-a[i]) - (3.0/h[i-1])*(a[i]-a[i-1]);
    }

    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    for(int i = 1; i < n; i++){
        l[i] = 2.0*(x[i+1] - x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
    }

    l[n] = 1.0;
    z[n] = 0.0;
    c[n] = 0.0;

    for(int j = n-1; j>= 0; j -= 1){
        c[j] = z[j] - mu[j]*c[j+1];
        b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3.0;
        d[j] = (c[j+1]-c[j])/(3*h[j]);
    }
   
       
}

int main(){

    int baseline = 50;

    int n = 1328;

    int r = 0;

    double x[n+1];
    double a[n+1];

    ifstream source;
    source.open("testdata.dat");
    
    for(string line; getline(source, line); r++){
        istringstream in(line);
        double xpt, ypt;
        in >> xpt;
        in >> ypt;
        //cout << xpt << "\t" << ypt << endl;
        x[r] = xpt;
        a[r] = ypt;
    }

    // subtract baseline
    for(int i = 0; i < 1329; i++)
        a[i] -= baseline;


    double h[n+1];
    double alpha[n+1];

    double l[n+1];
    double mu[n+1];
    double z[n+1];
    double b[n+1];
    double c[n+1];
    double d[n+1];    

    getCoefficients(n, x, a, b, c, d);

    for(int j = 0; j < n; j++){
        //cout << j << ": " << a[j] << "\t" << b[j] << "\t" << c[j] << "\t" << d[j] << endl;
    }

    double currX = x[n];// x[0];
    int partIndex;    

    
    sg_filter(17, 5, x, a);

    // generate points
    for(int i = 1; i <= 10000; i++){
        
        for(int j = 0; j < n; j++){
            if(currX <= x[j] && currX >= x[j+1]){
                partIndex = j;
                
                //cout << "partIndex: " << partIndex << endl;
                break;
            }       
        }

        double total_y = a[partIndex] + b[partIndex]*(currX-x[partIndex]) + c[partIndex]*pow((currX-x[partIndex]),2.0) + d[partIndex]*pow((currX-x[partIndex]), 3.0);

        cout << currX << "\t" << total_y << endl;
                

        currX = currX - (x[n]-x[0])/10000.0;
        
    }

    get_peaks(x, a);

    //boxcar_filter(5, x, a);
 
    //sg_filter(17, 5, x, a);   
}
