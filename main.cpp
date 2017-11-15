#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

double f_x(double currX, double x[1329], double a[1329], double b[1329], double c[1329], double d[1329]){

    int partIndex = 0;

    for(int j = 0; j < 1329-1; j++){
        if(currX <= x[j] && currX >= x[j+1]){
            partIndex = j;
            break;
        }
    }

    double total_y = a[partIndex] + b[partIndex]*(currX-x[partIndex]) + c[partIndex]*pow((currX-x[partIndex]),2.0) + d[partIndex]*pow((currX-x[partIndex]), 3.0);

    return total_y;
}

/*double composite_trapezoid(double a, double b, int n){
    double h = (b-a)/n;

    double XI0 = f_x(a) + f_x(b);
    double XI1 = 0;
    double XI2 = 0;

    for(int i = 1; i < n; i++){
        double X = a + i*h;
        cout << "X: " << X << endl;
        cout << "f(X): " << f_x(X) << endl;
        XI2 = XI2 + f_x(X);
    }

    XI1 = h*(XI0+2.0*XI2)/2.0;

    return XI1;
}


double composite_simpsons(double a, double b, int n){
    double h = (b-a)/n;

    double XI0 = f_x(a) + f_x(b);
    double XI1 = 0;
    double XI2 = 0;

    for(int i = 1; i < n; i++){
        double X = a + i*h;
        if(i%2 == 0)
            XI2 = XI2 + f(X);
        else
            XI1 = XI1 + f(X);
    }

    XI1 = h*(XI0+2.0*XI2+4.0*XI1)/3.0;

    return XI1;
}

double romberg(double a, double b, int n){

    double R[3][n+1];
    double h = b-a;

    R[1][1]=(h/2)*(f(a) + f(b));
    cout << R[1][1] << endl;

    for(int i = 2; i <= n; i++){

        double approx = 0;

        for(int k = 1; k <= pow(2.0, i-2.0); k++){
            approx = approx + f_x(a + (k-0.5)*h);
        }

        R[2][1] = 0.5 * ( R[1][1] + h*approx);

        cout << R[2][1] << " ";

        for(int j = 2; j <= i; j++){
            R[2][j] = R[2][j-1] + (R[2][j-1] - R[1][j-1])/(pow(4.0, j-1)-1.0);
            cout << R[2][j] << " ";
        }
        cout << endl;

        h = h/2.0;

        for (int j = 1; j <= i; j++)
            R[1][j] = R[2][j];
    }

    return R[2][n];

}

double adaptiveQuadrature(double a, double b, double TOL, int N){
    double APP = 0.0;
    int i = 1;
    double TOL_arr[N];
    double a_arr[N];
    double h_arr[N];
    double FA_arr[N];
    double FC_arr[N];
    double FB_arr[N];
    double S_arr[N];
    double L_arr[N];


    TOL_arr[i] = 10.0*TOL;
    a_arr[i] = a;
    h_arr[i] = (b-a)/2.0;
    FA_arr[i] = f_x(a);
    FC_arr[i] = f_x(a+h_arr[i]);
    FB_arr[i] = f_x(b);
    S_arr[i] = h_arr[i]*(FA_arr[i] + 4*FC_arr[i] + FB_arr[i])/3.0;
    L_arr[i] = 1.0;

    while(i>0){

        double FD = f(a_arr[i] + h_arr[i]/2.0);
        double FE = f(a_arr[i] + 3*h_arr[i]/2.0);
        double S1 = h_arr[i]*(FA_arr[i]+ 4*FD + FC_arr[i])/6.0;

        double S2 = h_arr[i]*(FC_arr[i] + 4*FE+FB_arr[i])/6.0;
        double v1 = a_arr[i];
        double v2 = FA_arr[i];
        double v3 = FC_arr[i];
        double v4 = FB_arr[i];
        double v5 = h_arr[i];
        double v6 = TOL_arr[i];
        double v7 = S_arr[i];
        double v8 = L_arr[i];

        i = i - 1;

        if(abs(S1 + S2 - v7) < v6){

            APP = APP + (S1 + S2);

        }
        else{
            if(v8 >= N){
                cout << "level exceeded" << endl;
                return NULL;
            }
            else{
                i = i + 1;
                a_arr[i] = v1 + v5;
                FA_arr[i] = v3;
                FC_arr[i] = FE;
                FB_arr[i] = v4;
                h_arr[i] = v5/2.0;
                TOL_arr[i] = v6/2.0;
                S_arr[i] = S2;
                L_arr[i] = v8 + 1.0;

                i = i + 1;
                a_arr[i] = v1;
                FA_arr[i] = v2;
                FC_arr[i] = FD;
                FB_arr[i] = v3;
                h_arr[i] = h_arr[i-1];
                TOL_arr[i] = TOL_arr[i-1];
                S_arr[i] = S1;
                L_arr[i] = L_arr[i-1];
            }
        }
    }

    return APP;
}
*/

double bisection(double a, double b, double x[1329], double aC[1329], double bC[1329], double cC[1329], double dC[1329], double tol){

    cout << "Bisection Algorithm" << endl;

    int n = 1;

    double fa = f_x(a, x, aC, bC, cC, dC);

    cout << "Iteration\t" << "Value of x\t" << "Value of f(x)\t" << "Absolute Error\t" << "Relative Error" <<  endl;

    while( n <= 100 ){
        double p = a + (b-a)/2.0;

        cout << n << "\t";

        printf("%15f\t", p);

        double fp = f_x(p, x, aC, bC, cC, dC);

        //double abs_err = abs(truth - p);

        //printf("%15f\t%15f\t%15f\n", fp, abs_err, abs_err/abs(truth));

        if(fp == 0 || (b-a)/2.0 < tol){
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

void get_peaks(double x[1329], double a[1329], double b[1329], double c[1329], double d[1329], double tol, vector<vector<double>> &peaks, double shift){

    bool tms_set = false;

    double y_prev = a[0];

    for(int i = 1; i < 1329; i++){
        if(y_prev < 0 and a[i] > 0){
            cout << "start between " << x[i-1] << " and " << x[i] << endl;
            double pntA = bisection(x[i-1], x[i], x, a, b, c, d, tol);
            vector<double> temp;
            temp.push_back(pntA);
            peaks.push_back(temp);
            cout << "point A: " << pntA << endl;
        }
        if(y_prev > 0 and a[i] < 0){
            cout << "end between " << x[i-1] << " and " << x[i] << endl;
            double pntB = bisection(x[i-1], x[i], x, a, b, c, d, tol);

            if(! tms_set){
                shift = (pntB+peaks.back().front())/2.0;
                for(int j = 0; j < 1329; j++){
                    x[j] -= shift;
                }
                pntB -= shift;
                peaks.back().front() -= shift;
                tms_set = true;
            }

            peaks.back().push_back(pntB);
            cout << "point B: " << pntB << endl;
            peaks.back().push_back((pntB+peaks.back().front())/2.0);

        }
        y_prev = a[i];
    }

}


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

    double baseline = 1400.0;

    double tol = 1e-8;

    int n = 1328;

    int r = 0;

    double shift;

    vector<vector<double>> peaks;


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
                break;
            }       
        }

        double total_y = a[partIndex] + b[partIndex]*(currX-x[partIndex]) + c[partIndex]*pow((currX-x[partIndex]),2.0) + d[partIndex]*pow((currX-x[partIndex]), 3.0);

        cout << currX << "\t" << total_y << endl;
                

        currX = currX - (x[n]-x[0])/10000.0;
        
    }

    get_peaks(x, a, b, c, d, tol, peaks, shift);

    //boxcar_filter(5, x, a);
 
    //sg_filter(17, 5, x, a);   

    for(int i = 0; i < peaks.size(); i++){
        for(int j = 0; j < peaks.at(i).size(); j++)
            cout << peaks.at(i).at(j) << "\t";;

        cout << endl;
    }

    cout << "tms shift " << shift << endl;
}
