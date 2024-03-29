#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

//Function to calculate Poisson probability
double poisson(double mu, int k) {
    double poisson_value = (pow (mu, k) * exp(-mu)) / tgamma(k + 1);
    return poisson_value;
}

int main() {
    using namespace std;
    int N = 234;                //numbers in datensumme.txt
    double mu = 3.11538;        //mean of datensumme.txt
    vector<int> zaehler(11);    //counting vector

    //read numbers from datensumme.txt and count how often the values k appear
    ifstream fin("datensumme.txt");
    int n_i;
    for(int i = 0 ; i < N ; ++i) {
        fin >> n_i;
        zaehler[n_i] += 1;
    }
    fin.close();

    //write count and expectation from the Poisson distribution to hist.txt
    ofstream fout("histpoi.txt");
    for(unsigned int k = 0; k < zaehler.size(); ++k){
      cout << k << ": " << zaehler[k] << endl;
      fout << k << " " << zaehler[k] << " " 
      << N * poisson(mu, k) << endl;

    }
    fout.close();
}