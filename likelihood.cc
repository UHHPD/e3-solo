#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

//Function to calculate Poisson probability
double poisson(double mu, int k) {
    double poisson_value = (pow (mu, k) * exp(-mu)) / tgamma(k + 1);
    return poisson_value;
}

//Funtcion to calculate the likelihood
double prob(std::vector<int> daten, double mu){
  double likelihood = 1;
  for(int k : daten){
    likelihood *= poisson(mu, k);
  }
  return likelihood;
}

//Function to calculate the likelihood(k_i)
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
    double mu = 3.11538;                      //mean of datensumme.txt
    double likelihood_max = 0;                //random start value
    double mu_max = 1;                        //random start value
    double log_likelihood_delta_min = 1e6;    //random start value
    double mu_min = 1;                        //random start value
    double lambda;                            //likelihood ratio
    double log_lambda;                        //-2*log(lambda)
    double n_dof;                             //degrees of freedom/mean of Chi^2 square distribution
    double standard_dev_chi_sq;              //standard deviation 
    double rel_dev_chi_sq;                   //relative deviation

    //read numbers from datensumme.txt and save them a vector
    ifstream fin("datensumme.txt");
    int n_i;
    for(int i = 0 ; i < 234 ; ++i) {
        fin >> n_i;
        daten.push_back(n_i);
    }
    fin.close();

    //write L(mu) - mu from 0 to 6 - to likelihood.txt
    //find Maximum of L(mu)
    //write -2*log(L(mu)) - mu from 0 to 6 - to nll.txt
    //write -2*log(L(mu))-(-2*log(3.11538))
    //find minimum of this(^) function 
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

      if (likelihood_value > likelihood_max){
        mu_max = j;
        likelihood_max = likelihood_value;
      }

      if (log_likelihood_delta_value < log_likelihood_delta_min){
        mu_min = j;
        log_likelihood_delta_min = log_likelihood_delta_value;
      }
    }
    fout_like.close();
    fout_log_like.close();
    fout_log_like_delta.close();

    //Calculation of lambda
    lambda = prob(daten,mu) / prob_k (daten);
    log_lambda = -2 * log(lambda);

    //Mean, standard deviation, relative deviation of Chi^2 distribution
    n_dof = 233;
    standard_dev_chi_sq = sqrt(2*n_dof);
    rel_dev_chi_sq = (log_lambda - n_dof) / standard_dev_chi_sq;

    //cout/print results
    cout << "L(mu = 3.11538)= " << prob(daten, mu) << endl;
    cout << "Maximum of likelihood:" << endl;
    cout << "mu_max = " << mu_max << endl;
    cout << "maximum value = " << likelihood_max << endl;
    cout << "Minimum of deltanll.txt:" << endl;
    cout << "mu_min = " << mu_min << endl;
    cout << "minimum value = " << log_likelihood_delta_min << endl;
    cout << "Chi^2:" << endl;
    cout << "-2 log (Lambda)= " << log_lambda << endl;
    cout << "Mean = " << n_dof << endl;
    cout << "Standard deviation= " << standard_dev_chi_sq << endl;
    cout << "Relative deviation= " << rel_dev_chi_sq << endl;
}