#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <time.h>
#include <complex>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

#define pi 3.14159265359

void gaussian(int n, vector<vector<complex<double>>> A, vector<complex<double>> b){

    for(int i = 0; i < n; i++){
        A[i].push_back(b[i]);
    }

    for(int i = 1; i <= n-1; i++){
        int p = -1;

        for(int j = 1; j <= n; j++){
            for(int k = 1; k <= n+1; k++){
                //cout << A[j-1][k-1] << "\t";
            }
            //cout << endl << endl;;
        }

        for(int j = i; j <= n; j++){

            //cout << A[j-1][i-1] << endl;

           if(A[j-1][i-1] != complex<double>(0)){

                if(p == -1)
                    p = j;
                else{
                    if(j<p)
                        p = j;
                }
                break;
            }
        }
        if(p == -1){
                cout << "No unique solution exists." << endl;
                return;
        }
        if(p != i){
            for(int k = 1; k <= n+1; k++){
                complex<double> temp = A[p-1][k-1];
                A[p-1][k-1] = A[i-1][k-1];
                A[i-1][k-1] = temp;
            }
        }

        for (int j = i + 1; j <= n; j++){
            complex<double> m = A[j-1][i-1]/A[i-1][i-1];
            for(int k = 1; k <= n+1; k++){
                A[j-1][k-1] = A[j-1][k-1] - m*A[i-1][k-1];
            }
        }
        /*
        for(int j = 1; j <= n; j++){
            for(int k = 1; k <= n+1; k++){
                cout << A[j-1][k-1] << "\t";
            }
            cout << endl;
        }
        cout << endl;*/
    }

    if(A[n-1][n-1] == 0.0){
        cout << "No unique solution exists" << endl;
        return;
    }

    complex<double> x[n];



    x[n-1] = A[n-1][n]/A[n-1][n-1];

    for(int i = n-1; i >= 1; i--){

        complex<double> sum = 0.0;
        for (int j = i+1; j <= n; j++){
            sum += A[i-1][j-1]*x[j-1];
        }

        x[i-1] = (A[i-1][n+1-1] - sum)/A[i-1][i-1];
    }

    cout << endl;

    for (int i = 0; i < n; i++){
        cout << x[i] << "\n";
    }
    cout << endl;



    for(int i = 0; i < n; i++)
        A[i].pop_back();
}

void dftnmr(int nFilterPoints, vector<double> x, vector<double> a){

    vector<vector<complex<double>>> z(nFilterPoints, vector<complex<double>>(nFilterPoints));

    cout << "z allocated" << endl;

    for(int j = 0; j < nFilterPoints; j++)
        for(int k = 0; k < nFilterPoints; k++)
            z[j][k] = (pow(exp(-sqrt(complex<double>(-1, 0))*complex<double>(2*pi/nFilterPoints)), j*k))/sqrt(nFilterPoints);

    cout << "z filled" << endl;

    for(int i = 0; i < 4; i++){
        cout << i << endl;
        for(int j = 0; j < 4; j++)
            cout << z[i][j] << endl;
    }

    vector<complex<double>> c;
    for(int i = 0; i < nFilterPoints; i++){
        complex<double> sum = 0.0;
        for(int j = 0; j < nFilterPoints; j++){
            sum += z[i][j]*complex<double>(a[i]);
        }
        c.push_back(sum);
    }

    for(int i = 0; i < nFilterPoints; i++)
        cout << c[i] << endl;


    vector<vector<double>> g(nFilterPoints, vector<double>(nFilterPoints));

    for(int i = 0; i < nFilterPoints; i++)
        for(int j = 0; j < nFilterPoints; j++){
            complex<double> d = 0.0;
            if(i == j)
                d = complex<double>(1);
            g[i][j] = (j-i)*exp((-4.0*log(2.0)*i*j)/(pow(nFilterPoints, 3.0/2.0)));
        }
    for(int i = 0; i < nFilterPoints; i++){
        complex<double> sum = 0.0;
        for(int j = 0; j < nFilterPoints; j++){
            sum += g[i][j]*(c[j]);
        }
        c[i] = (sum);
        cout << c[i] << endl;
    }

    gaussian(nFilterPoints, z, c);

}

/*
void dftnmr_inv2(int nFilterPoints, vector<double> x, vector<double> a){

    gsl_matrix_complex *z = gsl_matrix_alloc(nFilterPoints, nFilterPoints);

    gsl_matrix_complex *g = gsl_matrix_alloc(nFilterPoints, nFilterPoints);

    gsl_matrix_complex *b = gsl_matrix_alloc(nFilterPoints, nFilterPoints);

    gsl_permutation_complex *p = gsl_permutation_alloc (nFilterPoints);

    gsl_vector_complex *y = gsl_vector_alloc(nFilterPoints);

    gsl_vector_complex *c = gsl_vector_alloc(nFilterPoints);

    int i, j;

    for(int i = 0; i < nFilterPoints; i++){

        gsl_matrix_complex_set( y, i, j, a[i]);

    }


    for ( i=0; i<nFilterPoints; i++ )
        for ( j=0; j<nFilterPoints; j++ ){
            gsl_complex d1 = gsl_complex_rect(0.0, 1.0):
            gsl_complex a1 = gsl_complex_rect(-2*pi/nFilterPoints, 0.0):
            gsl_complex b1 = gsl_complex_rect(j*k, 0.0):
            gsl_complex c1 = gsl_complex_rect(sqrt(nFilterPoints), 0.0):

            gsl_complex setTo = gsl_complex_div(gsl_complex_pow(gsl_complex_exp(gsl_complex_mult(a1*d1)), b1), c1);

            gsl_matrix_complex_set( z, i, j, setTo);//(double) 1 / ( (double) (i+j-1)) );
        }
    for ( i=0; i<nFilterPoints; i++ )
        for ( j=0; j<nFilterPoints; j++ ){

            gsl_complex a1 = gsl_complex_rect(0.0, 0.0);

            gsl_complex b1 = gsl_complex_rect(exp((-4.0*log(2.0)*i*j)/(pow(nFilterPoints, 3.0/2.0))), 0.0);
            if(i == j) 
                a1 = gsl_complex_rect(1.0, 0.0);           
        }

    gsl_complex sum = gsl_complex_rect(0.0, 0.0);

    for ( i=0; i<nFilterPoints; i++ ){
        for ( j=0; j<nFilterPoints; j++ ){
            sum = gsl_complex_sum(sum, gsl_complex_mult(gsl_matrix_complex_get(z, i, j), gsl_vector_complex_get(c, i)));
        }
        gsl_vector_complex_set(c, i, sum);
    }


    sum = gsl_complex_rect(0.0, 0.0);
    for ( i=0; i<nFilterPoints; i++ ){
        for ( j=0; j<nFilterPoints; j++ ){
            sum = gsl_complex_sum(sum, gsl_complex_mult(gsl_matrix_complex_get(g, i, j), gsl_vector_complex_get(c, i)));
        }
        gsl_vector_complex_set(c, i, sum);
    }



    gsl_linalg_complex_LU_decomp (a, p, &permutation_sign);

    gsl_linalg_complex_LU_invert (a, p, b);

    gsl_matrix_complex_free(z);
    gsl_matrix_complex_free(b);
    gsl_permutation_complex_free (p);

}
*/
void dftnmr_inv(int nFilterPoints, vector<double> x, vector<double> a){
    // compound_1.dat
    
    vector<vector<complex<double>>> z(nFilterPoints, vector<complex<double>>(nFilterPoints));



    //cout << "z allocated" << endl;

    for(int j = 0; j < nFilterPoints; j++)
        for(int k = 0; k < nFilterPoints; k++)
            z[j][k] = (pow(exp(-(complex<double>(0, 1))*complex<double>(2*pi/nFilterPoints)), j*k))/sqrt(nFilterPoints);
            
    //cout << "z filled" << endl;

    //for(int i = 0; i < 4; i++){
        //cout << i << endl;
        //for(int j = 0; j < 4; j++)
            //cout << z[i][j] << endl;
    //}   

    vector<complex<double>> c;
    for(int i = 0; i < nFilterPoints; i++){
        complex<double> sum = 0.0; 
        for(int j = 0; j < nFilterPoints; j++){
            sum += z[i][j]*complex<double>(a[j]);   
        }
        c.push_back(sum);
    }

    //for(int i = 0; i < nFilterPoints; i++)
        //cout << c[i] << endl;

    
    vector<vector<double>> g(nFilterPoints, vector<double>(nFilterPoints));

    double d;

    for(int i = 0; i < nFilterPoints; i++)
        for(int j = 0; j < nFilterPoints; j++){

            if(i == j){
                d = 1;
            }
            else{
                d = 0;
            }

            g[i][j] = (d)*exp((-4.0*log(2.0)*i*j)/(pow(nFilterPoints, 3.0/2.0)));
        }

    for(int i = 0; i < nFilterPoints; i++){
        complex<double> sum = 0.0;
        for(int j = 0; j < nFilterPoints; j++){
            sum += g[i][j]*(c[j]);
        }
        c[i] = (sum);
        //cout << c[i] << endl;
    }

    vector<vector<complex<double>>> zinv(nFilterPoints, vector<complex<double>>(nFilterPoints));

    for(int j = 0; j < nFilterPoints; j++)
        for(int k = 0; k < nFilterPoints; k++)
            zinv[j][k] = conj(z[j][k]);

    for(int i = 0; i < nFilterPoints; i++){
        complex<double> sum = 0.0;
        for(int j = 0; j < nFilterPoints; j++){
            sum += zinv[i][j]*(c[j]);
        }

        //cout << sum << endl;
        a[i] = (sum).real();
        cout << x[i] << "\t" << a[i] << endl;
    }

}


void sg_filter(int nFilterPoints, int nPasses, vector<double> &x, vector<double> &a){

    double c5[5] = {-3, 12, 17, 12, -3};
    double c11[11] = {-36, 9, 44, 69, 84, 89, 84, 69, 44, 9, -36};
    double c17[17] = {-21, -6, 7, 18, 27, 34, 39, 42, 43, 42, 39, 34, 27, 18, 7, -6, -21};

    int m = (nFilterPoints - 1)/2;


    for(int p = 1; p <= nPasses; p++)
        for(int j = m; j < a.size()-m; j++){
            double sumX = 0.0;
            double sumA = 0.0;
            int count = 0;
            int N;
            if(nFilterPoints == 5)
                N = 35.0;
            else if(nFilterPoints == 11)
                N = 429.0;
            else if(nFilterPoints == 17)
                N = 323.0;

            for(int i = -m; i <= m; i++){

                if(nFilterPoints == 5){
                    sumA += c5[m+i]*a[i+j];
                }

                else if(nFilterPoints == 11){
                    sumA += c11[m+i]*a[i+j];
                }

                else if(nFilterPoints == 17){
                    sumA += c17[m+i]*a[i+j];
                }
                count++;
            }
            a[j] = sumA/N;
    }

}

int main(){

    ifstream source;
    source.open("testdata.dat");

    vector<double> x;
    vector<double> a;

    for(string line; getline(source, line);){
        istringstream in(line);
        double xpt, ypt;
        in >> xpt;
        in >> ypt;
        x.push_back(xpt);
        a.push_back(ypt);

        //cout << xpt << endl;

    }

    //cout << "points pushed" << endl;

    dftnmr_inv(x.size(), x, a);   

}


