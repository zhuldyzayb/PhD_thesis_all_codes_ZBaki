//
//  main.cpp
//  MALA sampling
//
//  Created by Zhuldyzay Baki on 19/09/2023.
//

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include "random.hpp"
#include "vector.hpp"
#include "pattern.hpp"
#include "cox2d.hpp"
#include "sampler.hpp"


using namespace std;
vector<vector<std::string>> readData(string name){ // just to read csv files
vector<vector<string>> data;
vector<string> row;
string line, word;
fstream file (name, ios::in);
    if(file.is_open()){
        while (getline(file, line)){
            row.clear();
            stringstream str(line);
            while (getline(str, word, ',')){
                row.push_back(word);
                data.push_back(row);
            }
        }}

    else { cout<<"Could not open the file\n";}
return(data);
}


int main() {
    whran.setSeed(789, 654, 123); // just to set seed for randomization
    
    //  TRUE PARAMETERS AND DATA
    // define csv files
     string stIDs="/Users/zhuldyzaybaki/Desktop/PhD/Cpp_codes/MALA sampling/Data/dataIDs.csv"; // IDs for filtering within the gas field W_s\times W_T
     string gasProd="/Users/zhuldyzaybaki/Desktop/PhD/Cpp_codes/MALA sampling/Data/dataVst.csv"; // V(s,t)
     string pressVal="/Users/zhuldyzaybaki/Desktop/PhD/Cpp_codes/MALA sampling/Data/dataMst.csv"; // m(s,t)
    // vectorised earthquake data: first Neq elements are x-coord, second Neq elements y-coord and the las t-values
    // x and y coordinates correspond to UTM zone 31, t is 0 for 1995 and 26 for 2021
     string eqPat="/Users/zhuldyzaybaki/Desktop/PhD/Cpp_codes/MALA sampling/Data/dataEQs.csv";
    // vectorised (by row) Q matrix of size NS by NS (Ns = number of sp locs)
     string Qmat="/Users/zhuldyzaybaki/Desktop/PhD/Cpp_codes/MALA sampling/Data/dataQmat.csv";
     
     vector<vector<string>> file1 = readData(stIDs);
     vector<vector<string>> file2 = readData(gasProd);
     vector<vector<string>> file3 = readData(pressVal);
     vector<vector<string>> file4 = readData(Qmat);
     vector<vector<string>> file5 = readData(eqPat);
     
     int N = ( (int) file1.size())-1;
     // define parameters
     double eta = -6.248427;
     double theta2 = 5.340544;
     double alpha = 0.095872;
     double beta = 0.1120963;
     double gam0 = alpha/beta;
     double sigma2 = 51.4089;
     double tausq = 22.06731;
     double tau = sqrt(tausq);
     double mu = 19.72706;
     
     // define sp-temp space
     double minx, maxx, miny, maxy, mint, maxt;
     int nx, ny, nt;
    // this data is taken from xrange and yrange after pixelation of Groningen shape file (eps =2km)
     minx = 737.2550; maxx = 771.5003; // min and max w
     miny = 5889.884; maxy = 5931.234; // min and max h
     mint = 0.0; maxt = 27.0;
     nx = 18; ny = 22; nt = 27;
     int sz = nx * ny * nt;
     if(N!=sz) {printf("Error: data input size does not match the number of grid locations.");}
     double delta = 1.0;
     double theta1 = eta + alpha*alpha*sigma2;
    
    
    
    int i;
    // READ DATA INTO VECTORS
    int *cInd = new int[N];
    for( i=0; i<N; i++ ) { cInd[i] = stoi(file1[i+1][0]); }
    double *cprod = new double[N];
    for( i=0; i<N; i++ ) { cprod[i] = stod(file2[i+1][0]); }
    double *cpress = new double[N];
    for( i=0; i<N; i++ ) { cpress[i] = stod(file3[i+1][0]); }
    double *cQ = new double[ (int) file4.size()-1 ];
    for( i=0; i<((int) file4.size()-1); i++ ) { cQ[i] = stod(file4[i+1][0]); }
    int nEqs = ( (int) file5.size() - 1 )/3;
    double *cEqs = new double[ (int) file5.size() - 1 ];
    for( i=0; i<((int) file5.size() - 1); i++ ) { cEqs[i] = stod(file5[i+1][0]); }
    
    // DEFINE PARTS
    double *cbase = new double[N]; // base rate
    for( i=0; i<N; i++){ cbase[i] = exp( theta1 + theta2*cprod[i] ); }
    // SIMULATE Y(s,t)
    double factor = 1 - exp(-2*mu*delta);
    double *cY = new double[N];  // Y(s,t)
    for( i=0; i<N; i++){
        if(i<nx*ny) { cY[i] = whran.normal(0.0, 1.0); }
        else { cY[i] = whran.normal(0.0, pow(factor, 0.5)); }
    }
    // SIMULATE X(s,t)
    double *cnoise = new double[N];  // X(s,t)
    for( i=0; i<N; i++){ cnoise[i] = cpress[i] + whran.normal(0, pow(sigma2, 0.5)); }
    
    // CALCULATE -\log(\Gamma(s,t))
    double *cExp = new double[N]; // -log Gamma(s,t)
    for( i=0; i<N; i++) {
        int k = i/(nx*ny);
        if(k==0) { cExp[i] = -log(gam0); }
        else { cExp[i] = -log( exp(-cExp[i - nx*ny]) + alpha*delta) - alpha*(cnoise[i] - cnoise[i-nx*ny]); }
    }
    
    // DEFINE THE MODEL
    Cox *richter = new CoxRateState( alpha, sigma2, tau, mu, nx, ny, nt,  minx, miny, mint, maxx, maxy, maxt, cpress );
    richter->report();
    // DEFINE PROPOSAL VARIANCES
    double *hh = new double[2];
    hh[0] = 0.02;
    hh[1] = 0.01;
    // DEFINE SAMPLER
    Sampler *sampler = new Sampler( richter, hh);
    sampler->report();
    
    // CONVERT THE OBSERVED DATA INTO PATTERN
    Pattern D;
    int nD;
    for( i=0; i<nEqs; i++){ D.add( new MarkedEvent(cEqs[i], cEqs[i+nEqs], cEqs[i+2*nEqs])); }
    nD = D.getSize();
    cout<<"Number of observed earthquakes: "<<nD<<endl;
   
    // CREATE \Lambda(s,t) for all (s,t)
    richter->createLambda( D, cExp, cbase, cnoise, cY, cQ, cInd, gam0, tau, mu, N);
    
    // ADRESSEs FOR SAVING OUTPUTS
    string adress = "/Users/zhuldyzaybaki/Desktop/PhD/Cpp_codes/MALA sampling/Sampled/";
    string fname; // for sampled random fields
    string ftrack; // for trace plots
        
    
    ofstream file;
    // DEFINE THE NUMBER OF ITERATIONS
    int nn = 5000;
    double *aa = new double[1];
    aa[0] = 0.0;
    // VECTORS FOR MEAN AND ACCEPTANCE IDENTIFIER AT EACH ITERATION (for trace plots)
    double *accept = new double[nn];
    double *mE = new double[nn];
    double *mY = new double[nn];
    
    for(int i=0; i<nn; i++){
        sampler->MALA_Step(cQ, aa, tau, mu); // one step sample
        accept[i] = aa[0];
        aa[0] = 0.0;
        double meanE = 0.0;
        double meanY = 0.0;
        // WRITE THE RESULTING RANDOM FIELD INTO SPECIFIED DIRECTORY
        fname = adress + "sample" + to_string(i+1) +".txt";
        file.open(fname);
        file<<"E(s,t); Y(s,t); -log[Gamma(s,t)]; \n";
        for(int i=0; i<N; i++){
            int tPos = i/(nx*ny);
            int xPos = i%nx;
            int yPos = i/nx - tPos * ny;
            file << richter->getRandom(xPos, yPos, tPos) - cpress[i]<<"; "<<richter->getSpRandom(xPos, yPos, tPos)<<"; "<<richter->getExponent(xPos, yPos, tPos)<<"\n";
                
            meanE += richter->getRandom(xPos, yPos, tPos) - cpress[i];
            meanY += richter->getSpRandom(xPos, yPos, tPos);
        }
    file.close();
    mE[i] = meanE/sz;
    mY[i] = meanY/sz;
    }
    // WRITE MEANS DATA INTO SPECIFIED DIRECTORY (trace plots)
    double vv = 0;
    ftrack = adress + "tracking.txt";
    file.open(ftrack);
    file<<"Mean E(s,t); Mean Y(s,t); Acceptance \n";
    for(int i=0; i<nn; i++){
        vv += accept[i];
        file<<mE[i]<<"; "<<mY[i]<<"; "<<accept[i]<<"\n";
        }
    cout<<"The acceptance rate is: "<< vv/nn<<endl;
    
    
    return 0;
}
