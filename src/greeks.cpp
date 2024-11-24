#include "../hd/greeks.h"
#include "../hd/solver.h" // Inclut les solveurs Crank-Nicolson pour Européennes et Américaines
#include <armadillo> // Inclure Armadillo pour les calculs de vecteurs et matrices




// Constructeur qui prend un objet Solver en paramètre
Greeks::Greeks(Solver solver) : solver_(solver) {}

// Fonction pour calculer la dérivée première (utilisée pour Delta)
arma::vec Greeks::derivee(const arma::vec& source, double h) {
    arma::vec result = arma::zeros<arma::vec>(source.n_elem - 1);
    for (size_t i = 0; i < result.n_elem; ++i) {
        result(i) = (source(i + 1) - source(i)) / h;
    }
    return result;
}

// Fonction pour calculer la dérivée seconde (utilisée pour Gamma)
arma::vec Greeks::deriveeSeconde(const arma::vec& source, double h) {
    arma::vec result = arma::zeros<arma::vec>(source.n_elem - 2);
    for (size_t i = 0; i < result.n_elem; ++i) {
        result(i) = (source(i + 2) - 2 * source(i + 1) + source(i)) / (h * h);
    }
    return result;
}

// Delta pour Call Européen
double Greeks::Delta_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Call = solver_.CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas = 2 * S0 / dis_spac;

    arma::vec last_col(dis_spac + 1);
    for (size_t i = 0; i <= dis_spac; ++i) {
        last_col(i) = Call(i, dis_tmp - 1);
    }

    arma::vec delta_values = derivee(last_col, pas);
    return delta_values(dis_spac / 2);
}

// Delta pour Put Européen
double Greeks::Delta_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Put = solver_.CrankNicolson_Put_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas = 2 * S0 / dis_spac;

    arma::vec last_col(dis_spac + 1);
    for (size_t i = 0; i <= dis_spac; ++i) {
        last_col(i) = Put(i, dis_tmp - 1);
    }

    arma::vec delta_values = derivee(last_col, pas);
    return delta_values(dis_spac / 2);
}

// Delta pour Call Américain
double Greeks::Delta_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Call = solver_.CrankNicolson_Call_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas = 2 * S0 / dis_spac;

    arma::vec last_col(dis_spac + 1);
    for (size_t i = 0; i <= dis_spac; ++i) {
        last_col(i) = Call(i, dis_tmp - 1);
    }

    arma::vec delta_values = derivee(last_col, pas);
    return delta_values(dis_spac / 2);
}

// Delta pour Put Américain
double Greeks::Delta_Put_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Put = solver_.CrankNicolson_Put_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas = 2 * S0 / dis_spac;

    arma::vec last_col(dis_spac + 1);
    for (size_t i = 0; i <= dis_spac; ++i) {
        last_col(i) = Put(i, dis_tmp - 1);
    }

    arma::vec delta_values = derivee(last_col, pas);
    return delta_values(dis_spac / 2);
}

// Gamma pour Call Européen
double Greeks::Gamma_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Call = solver_.CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas = 2 * S0 / dis_spac;

    arma::vec last_col(dis_spac + 1);
    for (size_t i = 0; i <= dis_spac; ++i) {
        last_col(i) = Call(i, dis_tmp - 1);
    }

    arma::vec gamma_values = deriveeSeconde(last_col, pas);
    return gamma_values(dis_spac / 2);
}

// Gamma pour Put Européen
double Greeks::Gamma_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Put = solver_.CrankNicolson_Put_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas = 2 * S0 / dis_spac;

    arma::vec last_col(dis_spac + 1);
    for (size_t i = 0; i <= dis_spac; ++i) {
        last_col(i) = Put(i, dis_tmp - 1);
    }

    arma::vec gamma_values = deriveeSeconde(last_col, pas);
    return gamma_values(dis_spac / 2);
}

// Gamma pour Call Américain
double Greeks::Gamma_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Call = solver_.CrankNicolson_Call_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas = 2 * S0 / dis_spac;

    arma::vec last_col(dis_spac + 1);
    for (size_t i = 0; i <= dis_spac; ++i) {
        last_col(i) = Call(i, dis_tmp - 1);
    }

    arma::vec gamma_values = deriveeSeconde(last_col, pas);
    return gamma_values(dis_spac / 2);
}

// Gamma pour Put Américain
double Greeks::Gamma_Put_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Put = solver_.CrankNicolson_Put_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas = 2 * S0 / dis_spac;

    arma::vec last_col(dis_spac + 1);
    for (size_t i = 0; i <= dis_spac; ++i) {
        last_col(i) = Put(i, dis_tmp - 1);
    }

    arma::vec gamma_values = deriveeSeconde(last_col, pas);
    return gamma_values(dis_spac / 2);
}

// Theta pour Call Européen
double Greeks::Theta_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Call = solver_.CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas_t = T / dis_tmp;

    arma::vec time_row(dis_tmp + 1);
    for (size_t i = 0; i <= dis_tmp; ++i) {
        time_row(i) = Call(dis_spac / 2, i);
    }

    arma::vec theta_values = derivee(time_row, pas_t);
    return -theta_values(theta_values.n_elem - 1);
}

// Theta pour Put Européen
double Greeks::Theta_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Put = solver_.CrankNicolson_Put_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas_t = T / dis_tmp;

    arma::vec time_row(dis_tmp + 1);
    for (size_t i = 0; i <= dis_tmp; ++i) {
        time_row(i) = Put(dis_spac / 2, i);
    }

    arma::vec theta_values = derivee(time_row, pas_t);
    return -theta_values(theta_values.n_elem - 1);
}

// Theta pour Call Américain
double Greeks::Theta_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Call = solver_.CrankNicolson_Call_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas_t = T / dis_tmp;

    arma::vec time_row(dis_tmp + 1);
    for (size_t i = 0; i <= dis_tmp; ++i) {
        time_row(i) = Call(dis_spac / 2, i);
    }

    arma::vec theta_values = derivee(time_row, pas_t);
    return -theta_values(theta_values.n_elem - 1);
}

// Theta pour Put Américain
double Greeks::Theta_Put_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Put = solver_.CrankNicolson_Put_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas_t = T / dis_tmp;

    arma::vec time_row(dis_tmp + 1);
    for (size_t i = 0; i <= dis_tmp; ++i) {
        time_row(i) = Put(dis_spac / 2, i);
    }

    arma::vec theta_values = derivee(time_row, pas_t);
    return -theta_values(theta_values.n_elem - 1);
}

// Vega pour Call Européen
double Greeks::Vega_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Call = solver_.CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas_sigma = sigma / 100;

    // Calcul de Vega par perturbation de la volatilité
    double price_up = solver_.CrankNicolson_Call_European(S0, K, T, r, sigma + pas_sigma, dis_tmp, dis_spac)(dis_spac / 2, dis_tmp - 1);
    double price_down = solver_.CrankNicolson_Call_European(S0, K, T, r, sigma - pas_sigma, dis_tmp, dis_spac)(dis_spac / 2, dis_tmp - 1);

    return (price_up - price_down) / (2 * pas_sigma);
}

// Vega pour Put Européen
double Greeks::Vega_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Put = solver_.CrankNicolson_Put_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas_sigma = sigma / 100;

    // Calcul de Vega par perturbation de la volatilité
    double price_up = solver_.CrankNicolson_Put_European(S0, K, T, r, sigma + pas_sigma, dis_tmp, dis_spac)(dis_spac / 2, dis_tmp - 1);
    double price_down = solver_.CrankNicolson_Put_European(S0, K, T, r, sigma - pas_sigma, dis_tmp, dis_spac)(dis_spac / 2, dis_tmp - 1);

    return (price_up - price_down) / (2 * pas_sigma);
}

// Rho pour Call Européen
double Greeks::Rho_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Call = solver_.CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas_r = r / 100;

    // Calcul de Rho par perturbation du taux d'intérêt
    double price_up = solver_.CrankNicolson_Call_European(S0, K, T, r + pas_r, sigma, dis_tmp, dis_spac)(dis_spac / 2, dis_tmp - 1);
    double price_down = solver_.CrankNicolson_Call_European(S0, K, T, r - pas_r, sigma, dis_tmp, dis_spac)(dis_spac / 2, dis_tmp - 1);

    return (price_up - price_down) / (2 * pas_r);
}

// Rho pour Put Européen
double Greeks::Rho_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    auto Put = solver_.CrankNicolson_Put_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
    double pas_r = r / 100;

    // Calcul de Rho par perturbation du taux d'intérêt
    double price_up = solver_.CrankNicolson_Put_European(S0, K, T, r + pas_r, sigma, dis_tmp, dis_spac)(dis_spac / 2, dis_tmp - 1);
    double price_down = solver_.CrankNicolson_Put_European(S0, K, T, r - pas_r, sigma, dis_tmp, dis_spac)(dis_spac / 2, dis_tmp - 1);

    return (price_up - price_down) / (2 * pas_r);
}

// Black-Scholes Delta pour Call
double Greeks::Delta_BS_Call(double S0, double K, double T, double r, double sigma) {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return arma::normcdf(d1);
}

// Black-Scholes Delta pour Put
double Greeks::Delta_BS_Put(double S0, double K, double T, double r, double sigma) {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return arma::normcdf(d1) - 1;
}

// Black-Scholes Gamma
double Greeks::Gamma_BS(double S0, double K, double T, double r, double sigma) {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return arma::normpdf(d1) / (S0 * sigma * sqrt(T));
}

// Black-Scholes Theta pour Call
double Greeks::Theta_BS_Call(double S0, double K, double T, double r, double sigma) {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    double term1 = -S0 * arma::normpdf(d1) * sigma / (2 * sqrt(T));
    double term2 = r * K * exp(-r * T) * arma::normcdf(d2);
    return term1 - term2;
}

// Black-Scholes Theta pour Put
double Greeks::Theta_BS_Put(double S0, double K, double T, double r, double sigma) {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    double term1 = -S0 * arma::normpdf(d1) * sigma / (2 * sqrt(T));
    double term2 = r * K * exp(-r * T) * arma::normcdf(-d2);
    return term1 + term2;
}

// Black-Scholes Vega
double Greeks::Vega_BS(double S0, double K, double T, double r, double sigma) {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return S0 * arma::normpdf(d1) * sqrt(T);
}

// Black-Scholes Rho pour Call
double Greeks::Rho_BS_Call(double S0, double K, double T, double r, double sigma) {
    double d2 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T)) - sigma * sqrt(T);
    return K * T * exp(-r * T) * arma::normcdf(d2);
}

// Black-Scholes Rho pour Put
double Greeks::Rho_BS_Put(double S0, double K, double T, double r, double sigma) {
    double d2 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T)) - sigma * sqrt(T);
    return -K * T * exp(-r * T) * arma::normcdf(-d2);
}
