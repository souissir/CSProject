#pragma once

#include <armadillo> // Pour les calculs de vecteurs et matrices

class Greeks {
    public :

    Greeks(Solver solver);
    // Crank-Nicolson : Delta, Gamma, Theta, Vega, Rho pour Européennes et Américaines
    double Delta_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Delta_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Delta_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Delta_Put_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);

    double Gamma_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Gamma_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Gamma_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Gamma_Put_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);

    double Theta_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Theta_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Theta_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Theta_Put_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);

    double Vega_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Vega_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);

    double Rho_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Rho_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);

    // Black-Scholes pour comparaison (Européennes uniquement)
    double Delta_BS_Call(double S0, double K, double T, double r, double sigma);
    double Delta_BS_Put(double S0, double K, double T, double r, double sigma);

    double Gamma_BS(double S0, double K, double T, double r, double sigma);

    double Theta_BS_Call(double S0, double K, double T, double r, double sigma);
    double Theta_BS_Put(double S0, double K, double T, double r, double sigma);

    double Vega_BS(double S0, double K, double T, double r, double sigma);

    double Rho_BS_Call(double S0, double K, double T, double r, double sigma);
    double Rho_BS_Put(double S0, double K, double T, double r, double sigma);

    double Vega_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    double Rho_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    arma::vec derivee(const arma::vec& source, double h);
    arma::vec deriveeSeconde(const arma::vec& source, double h);


    private :
      Solver solver_;
};
