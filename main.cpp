#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <time.h>
using namespace std;

// prints results to output file 
void print_output(string output_file, double baseline, double tol, int filter, int filter_size, int filter_passes, int integration, string infile, double shift, vector<vector<double>>peaks, double time){
    ofstream out_file(output_file);
    if(out_file.is_open()){
        out_file << "\t\t\t-=> NMR ANALYSIS <=-\n" << endl;
        out_file << "Program Options" << endl;
        out_file << "=============================" << endl;
        out_file << "Baseline Adjustment\t:\t" << baseline << endl;
        out_file << "Tolerance\t:\t" << tol << endl;
        
        if(filter == 0)
            out_file << "No filtering" << endl;
        else if(filter == 2){
            out_file << "Savitzky-Golay Filtering" << endl;
            out_file << "SG Filter Size (Cyclic)\t:\t" << filter_size << endl;
            out_file << "SG Filter Passes\t\t:\t" << filter_passes << endl;
        }
        else if(filter == 1){
            out_file << "Boxcar Filtering" << endl;
            out_file << "Boxcar Size (Cyclic)\t:\t" << filter_size << endl;
            out_file << "Boxcar Passes\t\t:\t" << filter_passes << endl;
        }
        out_file << endl;

        out_file << "Integration Method" << endl;
        out_file << "===========================" << endl;
        if(integration == 0)
            out_file << "Simpson's Rule" << endl;
        else if(integration == 1)
            out_file << "Trapezoid Rule" << endl; 
        else if(integration == 2)
            out_file << "Romberg Integration" << endl;
        else if(integration == 3)
            out_file << "Adaptive Quadrature" << endl;
        else if(integration == 4)
            out_file << "Gaussian Quadrature" << endl;
        out_file << endl;        

        out_file << "Plot File Data" << endl;       
        out_file << "===========================" << endl;
        out_file << "File: " << infile << endl;
        out_file << "Plot shifted\t" << shift << " for the TMS calibration" << endl << endl << endl;

        out_file << "Peak   " << setw(10)<< "Begin" << setw(10) << "End" << setw(10) << "Location" << setw(10) << "Top" << setw(10) <<"Area" << setw(10) << "Hydrogens" << endl;
        out_file << setw(7) << "=======" << setw(10) << "========" << setw(10) << "========" << setw(10) << "========" << setw(10) << "========" << setw(10) << "========" << setw(10) << "========" << endl;  

        int size_peaks = peaks.size();
        
        int peak_count = 1;
        for(int i = peaks.size()-1; i >=0; i--){
            out_file  << setw(7) << peak_count;
            out_file << std::setw(10) << peaks.at(i).at(1);
            out_file << std::setw(10) << peaks.at(i).at(0);
            for(int j = 2; j < peaks.at(i).size(); j++)
                out_file << std::setw(10) << peaks.at(i).at(j);
            out_file << endl;
            peak_count++;
        }
        
        out_file << "Anaylsis took " << time << " seconds" << endl;   

        out_file.close();
    }
}

// read analysis specifications from nmr.in
void read_options(string &input_file, double &baseline, double &tol, int &filter, int &filter_size, int &filter_passes, int &integration, string &output_file, int &n_points){
    string line;
    ifstream options_file("nmr.in");
    if(options_file.is_open()){
        getline(options_file, input_file);
        
        getline(options_file, line);
        baseline = atof(line.c_str());

        getline(options_file, line);
        tol = atof(line.c_str());

        getline(options_file, line);
        filter = atoi(line.c_str());

        getline(options_file, line);
        filter_size = atoi(line.c_str());

        getline(options_file, line);
        filter_passes = atoi(line.c_str());

        getline(options_file, line);
        integration = atoi(line.c_str());
        
        getline(options_file, output_file);
   
        options_file.close();

        // count lines in input_file
        ifstream infile(input_file);
        int n_points = 0;       
        if(infile.is_open()){
            while(getline(infile, line)){
                n_points++;
            }
            infile.close();
            //cout << "n points " << n_points << endl;
        }
        else{
            cout << "input file not found" << endl;
        }
    }
    else{
        cout << "nmr.in not opened" << endl;
    }
    
}

// cubic spline function
double f_x(double currX, vector<double> &x, vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d){

    int partIndex = 0;

    for(int j = 0; j < x.size() - 1; j++){
        if(currX <= x.at(j) && currX >= x.at(j+1)){
            partIndex = j;
            break;
        }
    }

    double total_y = a.at(partIndex) + b.at(partIndex)*(currX-x.at(partIndex)) + c.at(partIndex)*pow(currX-x.at(partIndex), 2.0) + d.at(partIndex)*pow(currX-x.at(partIndex), 3.0);
    return total_y;
}

// derivative of cubic spline function
double df_x(double currX, vector<double> &x, vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d){

    int partIndex = 0;

    for(int j = 0; j < x.size() - 1; j++){
        if(currX <= x.at(j) && currX >= x.at(j+1)){
            partIndex = j;
            break;
        }
    }

    double total_y =  b.at(partIndex) + 2.0*c.at(partIndex)*pow(currX-x.at(partIndex), 1.0) + 3.0*d.at(partIndex)*pow(currX-x.at(partIndex), 2.0);
    return total_y;
}

// composite trapezoid method
double composite_trapezoid(double a, double b, int n, vector<double> &x, vector<double> &aC, vector<double> &bC, vector<double> &cC, vector<double> &dC){
    double h = (b-a)/n;

    double XI0 = f_x(a, x, aC, bC, cC, dC) + f_x(b, x, aC, bC, cC, dC);
    double XI1 = 0;
    double XI2 = 0;

    for(int i = 1; i < n; i++){
        double X = a + i*h;
        XI2 = XI2 + f_x(X, x, aC, bC, cC, dC);
    }

    XI1 = h*(XI0+2.0*XI2)/2.0;

    return XI1;
}

// composite Simpson's Method
double composite_simpsons(double a, double b, int n, vector<double> &x, vector<double> &aC, vector<double> &bC, vector<double> &cC, vector<double> &dC){

    double h = (b-a)/n;

    double XI0 = f_x(a, x, aC, bC, cC, dC) + f_x(b, x, aC, bC, cC, dC);
    double XI1 = 0;
    double XI2 = 0;

    for(int i = 1; i < n; i++){
        double X = a + i*h;
        if(i%2 == 0)
            XI2 = XI2 + f_x(X, x, aC, bC, cC, dC);
        else
            XI1 = XI1 + f_x(X, x, aC, bC, cC, dC);

        //cout << XI1 << endl;
    }

    XI1 = h*(XI0+2.0*XI2+4.0*XI1)/3.0;

    return XI1;
}

// Romberg Integration
double romberg(double a, double b, int n, vector<double> &x, vector<double> &aC, vector<double> &bC, vector<double> &cC, vector<double> &dC, double tol){

    double R[3][n+1];
    double h = b-a;

    R[1][1]=(h/2)*(f_x(a, x, aC, bC, cC, dC) + f_x(b,x, aC, bC, cC, dC));

    double prev = R[1][1];
    double diff = abs(0.0 - R[1][1]);

    for(int i = 2; i <= n && diff > tol; i++){

        double approx = 0;

        for(int k = 1; k <= pow(2.0, i-2.0); k++){
            approx = approx + f_x(a + (k-0.5)*h, x, aC, bC, cC, dC);
        }

        R[2][1] = 0.5 * ( R[1][1] + h*approx);

        for(int j = 2; j <= i; j++){
            R[2][j] = R[2][j-1] + (R[2][j-1] - R[1][j-1])/(pow(4.0, j-1)-1.0);
        }

        h = h/2.0;

        diff = abs(R[2][i]-prev);
        prev = R[2][i];

        for (int j = 1; j <= i; j++)
            R[1][j] = R[2][j];
    }

    return R[2][n];

}

// Adaptive Quadrature
double adaptiveQuadrature(double a, double b, double TOL, int N, vector<double> &x, vector<double> &aC, vector<double> &bC, vector<double> &cC, vector<double> &dC){

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
    FA_arr[i] = f_x(a, x, aC, bC, cC, dC);
    FC_arr[i] = f_x(a+h_arr[i], x, aC, bC, cC, dC);
    FB_arr[i] = f_x(b, x, aC, bC, cC, dC);
    S_arr[i] = h_arr[i]*(FA_arr[i] + 4*FC_arr[i] + FB_arr[i])/3.0;
    L_arr[i] = 1.0;

    while(i>0){

        double FD = f_x(a_arr[i] + h_arr[i]/2.0, x, aC, bC, cC, dC);
        double FE = f_x(a_arr[i] + 3*h_arr[i]/2.0, x, aC, bC, cC, dC);
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

            //cout << "APP "  << APP << endl;
        }
    }

    return APP;
}

// bistection method for finding intersection between splines and baseline
double bisection(double a, double b, vector<double> &x, vector<double> &aC, vector<double> &bC, vector<double> &cC, vector<double> &dC, double tol){

    int n = 1;

    double fa = f_x(a, x, aC, bC, cC, dC);

    while( n <= 100 ){
        double p = a + (b-a)/2.0;

        double fp = f_x(p, x, aC, bC, cC, dC);

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

// bisection algorithm for derivative of cubic splilne function
double dfxbisection(double a, double b, vector<double> &x, vector<double> &aC, vector<double> &bC, vector<double> &cC, vector<double> &dC, double tol){

    int n = 1;

    double fa = df_x(a, x, aC, bC, cC, dC);

    while( n <= 100 ){
        double p = a + (b-a)/2.0;

        double fp = df_x(p, x, aC, bC, cC, dC);

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

// shifts x points by maximum value of highest valued peak
void set_tms_shift(vector<double> &x, vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d, double tol, double &shift, double baseline){

    //subtract baseline
    for(int i = 0; i < x.size(); i++){
        a.at(i) -= baseline;
    }

    double y_prev = a.at(0);
    bool tms_set = false;
    double pntA = 0.0;

    int maxIndex;

    for(int i = 1; i < x.size() && !tms_set; i++){

        // get point a (where peak spline crosses baseline)
        if(y_prev < 0 and a.at(i) > 0){
            pntA = bisection(x.at(i-1), x.at(i), x, a, b, c, d, tol);
            maxIndex = i;
        }


        if(y_prev > 0 and a.at(i) < 0){//a[i] < 0){
            //get point B (where other end of peak crosses baseline)
            double pntB = bisection(x.at(i-1), x.at(i), x, a, b, c, d, tol);
            // set shift value and shift x values
            shift = x.at(maxIndex-1);
            for(int j = 0; j < a.size(); j++){
                x.at(j) -= shift;
            }
            tms_set = true;
        }
        // get maximum point in peak
        if(a.at(i) > y_prev)
            maxIndex = i;

        y_prev = a.at(i);   
    }

    //restore baseline
    for(int i = 0; i < x.size(); i++){
        a.at(i) += baseline;
    }

}

// gaussian integration when n = 5 
double gaussian(double pnta, double pntb, vector<double> &x, vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d){ 

    double r[5] = {0.9061798459, 0.5384693101, 0.0000, -0.5384693101, -0.9061798459};
    double coef[5] = {0.2369268885, 0.4786286705, 0.5688888889, 0.4786286705, 0.236926885};


    double sum = 0.0;
    for(int i = 0; i < 5; i++)
        sum += coef[i]*f_x(0.5*((pntb-pnta)*r[i]+pnta+pntb), x, a, b, c, d);

    return sum;
}

// gets peak data for each peak located
void get_peaks(vector<double> &x, vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d, double tol, vector<vector<double>> &peaks, double baseline, double integration){

    // subtract baseline
    for(int i = 0; i < a.size(); i++){
        a.at(i) -= baseline;
    }

    bool tms_set = false;

    double y_prev = a.at(0);

    double max_a = 0.0;

    for(int i = 1; i < a.size(); i++){
        if(a.at(i) > max_a){
            max_a = a.at(i);
        }

        if(y_prev < 0 and a.at(i) > 0){
            // create new peak and add end point
            double pntA = bisection(x[i-1], x[i], x, a, b, c, d, tol);
            vector<double> temp;
            temp.push_back(pntA);
            peaks.push_back(temp);
        }
        if(y_prev > 0 and a.at(i) < 0){

            // find start point and add to peak data
            double pntB = bisection(x.at(i-1), x.at(i), x, a, b, c, d, tol);
            peaks.back().push_back(pntB);
    
            // add midpoint of peak
            peaks.back().push_back((pntB+peaks.back().front())/2.0);
            // add maximum value in peak
            peaks.back().push_back(max_a);
      
            // get peak area
            double area;
            if(integration == 0)
                area = composite_simpsons(pntB, peaks.back().front(), 10, x, a, b, c, d);
            else if(integration == 1)
                area = composite_trapezoid(pntB, peaks.back().front(), 10, x, a, b, c, d);   
            else if(integration == 2)
                area = romberg(pntB, peaks.back().front(), 7, x, a, b, c, d, tol);
            else if(integration == 3)
                area = adaptiveQuadrature(pntB, peaks.back().front(), tol, 30, x, a, b, c, d);
            else if(integration == 4)
                area = gaussian(pntB, peaks.back().front(), x, a, b, c, d);

            peaks.back().push_back(abs(area));
            max_a = 0.0;
        }
        y_prev = a.at(i);
    }

    // get hydrogen estimates
    
    double min_area = peaks.at(0).back();
    for (int i = 1; i < peaks.size(); i++){
        if(peaks.at(i).back() < min_area){
            min_area = peaks.at(i).back();
        }
    }

    for (int i = 0; i < peaks.size(); i++){
        peaks.at(i).push_back(floor(peaks.at(i).back()/min_area));
    }


    // restore baseline
    for(int i = 0; i < a.size(); i++){
        a.at(i) += baseline;
    }
}

// SG Filter for n = 5, 11, 17
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

// boxcar filter for odd values of n
void boxcar_filter(int nFilterPoints, int nPasses, vector<double>&x, vector<double> &a){
    
   for(int m = 0; m < nPasses; m++)
    for(int i = (nFilterPoints-1)/2; i < a.size()-(nFilterPoints-1)/2; i++){
        double sumX = 0.0;
        double sumA = 0.0;
        for(int j = i - (nFilterPoints-1)/2; j <= i + (nFilterPoints-1)/2; j++){
            sumA += a.at(j);
        }
        a.at(i) = sumA/nFilterPoints;
    }
    
}

// calculates cubic spline coefficients
void getCoefficients(int n, vector<double>&x, vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d){

    n = x.size()-1;

    //cout << n << endl;

    vector<double> h;
    vector<double> alpha;
    vector<double> l;
    vector<double> mu;
    vector<double> z;

    for(int i = 0; i < n+1; i++){
        h.push_back(0.0);
        alpha.push_back(0.0);
        l.push_back(0.0);
        mu.push_back(0.0);
        z.push_back(0.0);
    }
    
    for(int i = 0; i < n; i++){
        h.at(i) = x.at(i+1)-x.at(i);
    }

    for(int i = 1; i < n; i++){
        alpha.at(i) = (3.0/h.at(i))*(a.at(i+1)-a.at(i))-(3.0/h.at(i-1))*(a.at(i)-a.at(i-1));
    }

    l.at(0) = 1.0;
    mu.at(0) = 0.0;
    z.at(0) = 0.0;

    for(int i = 1; i < n; i++){
        l.at(i)=2.0*(x.at(i+1)-x.at(i-1))-h.at(i-1)*mu.at(i-1);
        mu.at(i) = h.at(i)/l.at(i);
        z.at(i) = (alpha.at(i) - h.at(i-1)*z.at(i-1))/l.at(i);
    }

    l.at(n) = 1.0;
    z.at(n) = 0.0;
    c.at(n) = 0.0;

    for(int j = n-1; j>= 0; j -= 1){
        c.at(j) = z.at(j) - mu.at(j)*c.at(j+1);
        b.at(j) = (a.at(j+1)-a.at(j))/h.at(j) - h.at(j)*(c.at(j+1)+2*c.at(j))/3.0;
        d.at(j) = (c.at(j+1)-c.at(j))/(3*h.at(j));
    }
   
       
}

int main(){

    double baseline;
    double tol;
    int n;
    int r = 0;
    double shift = 0.0;
    string input_file, output_file;
    int filter, filter_size, filter_passes, integration;
    clock_t start = clock(), diff;
    
    // read options from nmr.in file
    read_options(input_file, baseline, tol, filter, filter_size, filter_passes, integration, output_file, n);

    n = n-1;

    vector<vector<double>> peaks;
    vector<double> x;
    vector<double> a;

    // read data from input file
    ifstream source;
    source.open(input_file);
    
    for(string line; getline(source, line); r++){
        istringstream in(line);
        double xpt, ypt;
        in >> xpt;
        in >> ypt;
        x.push_back(xpt);
        a.push_back(ypt);
    }
    
    vector<double> h;
    vector<double> alpha;
    vector<double> l;
    vector<double> mu;
    vector<double> z;
    vector<double> b;
    vector<double> c;
    vector<double> d;

    for(int i = 0; i < n+1; i++){

        b.push_back(1.0);
        c.push_back(1.0);
        d.push_back(1.0);
    }

    // shift points relative to tms peak
    set_tms_shift(x, a, b, c, d, tol, shift, baseline);

    // filter shifted points
    if(filter == 1)
        boxcar_filter(filter_size, filter_passes, x, a);//(5, 5, x, a);
    else if(filter == 2)
        sg_filter(filter_size, filter_passes, x, a);

    // get cubic spline
    getCoefficients(n, x, a, b, c, d);

    // get cubic spline peaks
    get_peaks(x, a, b, c, d, tol, peaks, baseline, integration);

    n = x.size() - 1;

    double currX = x.at(n);
    int partIndex;    

    // generate points for graphing
    for(int i = 1; i <= 10000; i++){
        
        for(int j = 0; j < n; j++){
            if(currX <= x.at(j) && currX >= x.at(j+1)){
                partIndex = j;
                break;
            }       
        }

        double total_y = a.at(partIndex) + b.at(partIndex)*(currX-x.at(partIndex)) + c.at(partIndex)*pow(currX-x.at(partIndex), 2.0) + d.at(partIndex)*pow(currX-x.at(partIndex), 3.0);

        //uncomment for graphing points
        //cout << currX << "\t" << total_y << endl;        

        currX = currX -(x.at(n)-x.at(0))/10000.0;
    }


    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;

    print_output(output_file, baseline, tol, filter, filter_size, filter_passes, integration, input_file, shift, peaks, msec/1000.0);

}
