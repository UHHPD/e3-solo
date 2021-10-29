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
  return likelihood;
}

double prob_k(std::vector<int> daten){
  double likelihood = 1;
  for(int k: daten){
    likelihood *= poisson(k, k);
  }
  return likelihood;
}

int main() {
    using namespace std;
    vector<int> daten;
    double mu = 3.11538;
    double log_likelihood_delta_min = 1e6;        //random start value
    double mu_min = 1;        //random start value
    vector<double> interval;

    ifstream fin("datensumme.txt");
    int n_i;
    for(int i = 0 ; i < 234 ; ++i) {
        fin >> n_i;
        daten.push_back(n_i);
    }
    fin.close();

    cout << "L(mu = 3.11538)= " << prob(daten, mu) << std::endl;

    ofstream fout_like("likelihood.txt");
    ofstream fout_log_like("nll.txt");
    ofstream fout_log_like_delta("deltanll.txt");
    for(double j = 0; j <= 6; j += 0.001){
      double likelihood_value = prob(daten, j);
      double log_likelihood_value = -2 * log(prob(daten, j));
      double log_likelihood_delta_value = log_likelihood_value - (-2 * log(prob(daten, mu)));

      fout_like << j << " " << likelihood_value << std::endl;
      fout_log_like << j << " " << log_likelihood_value << std::endl;
      fout_log_like_delta << j << " " << log_likelihood_delta_value << std::endl;
    }
    fout_like.close();
    fout_log_like.close();
    fout_log_like_delta.close();

    double lambda = prob(daten,mu) / prob_k (daten);
    double ratio_value = -2 * log(lambda);

    cout << "-2 log (Lambda)= " << ratio_value << std::endl;

    double n_dof = 233;
    double standard_dev_chi_sq = sqrt(2*n_dof);
    double rel_dev_chi_sq = (ratio_value - n_dof) / standard_dev_chi_sq;

    cout << "Chi^2:" << std::endl;
    cout << "Mean = " << n_dof << std::endl;
    cout << "Standard deviation= " << standard_dev_chi_sq << std::endl;
    cout << "Relative deviation= " << rel_dev_chi_sq << std::endl;
}