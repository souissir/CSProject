#include "greeks.h"
#include "solver.h" // Inclut les solveurs Crank-Nicolson pour Européennes et Américaines
#include <vector>
#include <cmath> // Pour log, sqrt, exp

namespace greeks {

    // Fonction de dérivée première
    std::vector<double> Derivee(const std::vector<double>& source, double h) {
        std::vector<double> result(source.size() - 1);
        for (size_t i = 0; i < result.size(); ++i) {
            result[i] = (source[i + 1] - source[i]) / h;
        }
        return result;
    }

    // Fonction de dérivée seconde
    std::vector<double> DeriveeSeconde(const std::vector<double>& source, double h) {
        std::vector<double> result(source.size() - 2);
        for (size_t i = 0; i < result.size(); ++i) {
            result[i] = (source[i + 2] - 2 * source[i + 1] + source[i]) / (h * h);
        }
        return result;
    }

    // Delta pour Call Européen
    double Delta_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        auto Call = CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
        double pas = 2 * S0 / dis_spac;

        std::vector<double> last_col(dis_spac + 1);
        for (size_t i = 0; i <= dis_spac; ++i) {
            last_col[i] = Call[i][dis_tmp - 1];
        }

        std::vector<double> delta_values = Derivee(last_col, pas);
        return delta_values[dis_spac / 2];
    }

    // Gamma pour Call Européen
    double Gamma_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        auto Call = CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
        double pas = 2 * S0 / dis_spac;

        std::vector<double> last_col(dis_spac + 1);
        for (size_t i = 0; i <= dis_spac; ++i) {
            last_col[i] = Call[i][dis_tmp - 1];
        }

        std::vector<double> gamma_values = DeriveeSeconde(last_col, pas);
        return gamma_values[dis_spac / 2];
    }

    // Theta pour Call Européen
    double Theta_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        auto Call = CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac);
        double pas_t = T / dis_tmp;

        std::vector<double> time_row(dis_tmp + 1);
        for (size_t i = 0; i <= dis_tmp; ++i) {
            time_row[i] = Call[dis_spac / 2][i];
        }

        std::vector<double> theta_values = Derivee(time_row, pas_t);
        return -theta_values[theta_values.size() - 1];
    }

    // Vega pour Call Européen (par variation de sigma)
    double Vega_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        double price_1 = CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac)[dis_spac / 2][dis_tmp - 1];
        double price_2 = CrankNicolson_Call_European(S0, K, T, r, sigma + 0.01, dis_tmp, dis_spac)[dis_spac / 2][dis_tmp - 1];

        return (price_2 - price_1) / 0.01;
    }

    // Rho pour Call Européen (par variation de r)
    double Rho_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        double price_1 = CrankNicolson_Call_European(S0, K, T, r, sigma, dis_tmp, dis_spac)[dis_spac / 2][dis_tmp - 1];
        double price_2 = CrankNicolson_Call_European(S0, K, T, r + 0.01, sigma, dis_tmp, dis_spac)[dis_spac / 2][dis_tmp - 1];

        return (price_2 - price_1) / 0.01;
    }



    // Delta pour Call Américain
    double Delta_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        auto Call = CrankNicolson_Call_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
        double pas = 2 * S0 / dis_spac;

        std::vector<double> last_col(dis_spac + 1);
        for (size_t i = 0; i <= dis_spac; ++i) {
            last_col[i] = Call[i][dis_tmp - 1];
        }

        std::vector<double> delta_values = Derivee(last_col, pas);
        return delta_values[dis_spac / 2];
    }

    // Gamma pour Call Américain
    double Gamma_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        auto Call = CrankNicolson_Call_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
        double pas = 2 * S0 / dis_spac;

        std::vector<double> last_col(dis_spac + 1);
        for (size_t i = 0; i <= dis_spac; ++i) {
            last_col[i] = Call[i][dis_tmp - 1];
        }

        std::vector<double> gamma_values = DeriveeSeconde(last_col, pas);
        return gamma_values[dis_spac / 2];
    }

    // Theta pour Call Américain
    double Theta_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        auto Call = CrankNicolson_Call_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
        double pas_t = T / dis_tmp;

        std::vector<double> time_row(dis_tmp + 1);
        for (size_t i = 0; i <= dis_tmp; ++i) {
            time_row[i] = Call[dis_spac / 2][i];
        }

        std::vector<double> theta_values = Derivee(time_row, pas_t);
        return -theta_values[theta_values.size() - 1];
    }

    // Vega pour Call Américain
    double Vega_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        double price_1 = CrankNicolson_Call_American(S0, K, T, r, sigma, dis_tmp, dis_spac)[dis_spac / 2][dis_tmp - 1];
        double price_2 = CrankNicolson_Call_American(S0, K, T, r, sigma + 0.01, dis_tmp, dis_spac)[dis_spac / 2][dis_tmp - 1];

        return (price_2 - price_1) / 0.01;
    }

    // Rho pour Call Américain
    double Rho_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        double price_1 = CrankNicolson_Call_American(S0, K, T, r, sigma, dis_tmp, dis_spac)[dis_spac / 2][dis_tmp - 1];
        double price_2 = CrankNicolson_Call_American(S0, K, T, r + 0.01, sigma, dis_tmp, dis_spac)[dis_spac / 2][dis_tmp - 1];

        return (price_2 - price_1) / 0.01;
    }

    // Répéter pour Put Américain en remplaçant Call par Put
    double Delta_Put_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
        auto Put = CrankNicolson_Put_American(S0, K, T, r, sigma, dis_tmp, dis_spac);
        double pas = 2 * S0 / dis_spac;

        std::vector<double> last_col(dis_spac + 1);
        for (size_t i = 0; i <= dis_spac; ++i) {
            last_col[i] = Put[i][dis_tmp - 1];
        }

        std::vector<double> delta_values = Derivee(last_col, pas);
        return delta_values[dis_spac / 2];
    }
}

