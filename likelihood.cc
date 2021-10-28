#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

double poisson(double mu, int k) {
    double poisson_value = (pow (mu, k) * exp(-mu)) / tgamma(k + 1);
    return poisson_value;
}

double prob(std::vector<int> daten, double mu){
  double likelihood = 1;
  for(int k : daten){
    likelihood *= poisson(mu, k);
  }
}

int main() {
    using namespace std;
    vector<int> daten;
    double mu = 3.11538;
    double log_likelihood_delta_min = 1e6;        //random start value
    double mu_min = 1;        //random start value
    vector<double> interval;
    double mu_value;
    double delta_value;

    ifstream fin("datensumme.txt");
    int n_i;
    for(int i = 0 ; i < 234 ; ++i) {
        fin >> n_i;
        daten.push_back(n_i);
    }
    fin.close();

    cout << "L(mu = 3.11538): " << prob(daten, mu) << std::endl;

    ofstream fout_like("likelihood.txt");
    ofstream fout_log_like("nll.txt");
    ofstream fout_log_like_delta("deltanll.txt");
    for(double j = 0; j <= 6; j += 0.001){
      double likelihood_value = prob(daten, j);
      double log_likelihood_value = -2 * log(prob(daten, j));
      double log_likelihood_delta_value = log_likelihood_value - (-2 * log(prob(daten, mu)));

      if (log_likelihood_delta_value < log_likelihood_delta_min){
        log_likelihood_delta_min = log_likelihood_delta_value;
        mu_min = j;
      }

      fout_like << j << " " << likelihood_value << std::endl;
      fout_log_like << j << " " << log_likelihood_value << std::endl;
      fout_log_like_delta << j << " " << log_likelihood_delta_value << std::endl;
    }
    fout_like.close();
    fout_log_like.close();
    fout_log_like_delta.close();

    cout << "Mu min: " << mu_min << std::endl;
    cout << "Minimum delta log likelihood:  " << 
    log_likelihood_delta_min << std::endl;
  
}