#include "../hd/solver.h"
#include <algorithm>
#include <iostream> // For debugging or logging
#include <armadillo>

// Constructor
Solver::Solver(const std::string& optionType, const std::string& exerciseType,
               double maturity, double strike, double computationDate,
               const std::pair<double, double>& timeMeshParams,
               const std::pair<double, double>& spotMeshParams,
               double spotPrice, const std::vector<std::pair<double, double>>& riskFreeRate,
               double volatility)
    : optionType(optionType), exerciseType(exerciseType), maturity(maturity),
      strike(strike), computationDate(computationDate), spotPrice(spotPrice),
      timeMeshParams(timeMeshParams), spotMeshParams(spotMeshParams),
      riskFreeRate(riskFreeRate), volatility(volatility) {}

// Solve method
void Solver::solve() {
    setupGrid();
    applyBoundaryConditions();

    // Backward in time computation using the Crank-Nicholson method
    for (size_t t = timeSteps.n_elem - 1; t > 0; --t) {
        crankNicholsonStep(t);
    }
}

// Set up the grid
void Solver::setupGrid() {
    double S_min = strike / 3.0;       // Minimum spot price
    double S_max = strike * 2.0;       // Maximum spot price
    double dS = (S_max - S_min) / spotMeshParams.first;
    double dt = maturity / timeMeshParams.first;

    // Initialize asset prices and time steps using Armadillo
    assetPrices = arma::linspace(S_min, S_max, spotMeshParams.first + 1);
    timeSteps = arma::linspace(computationDate, maturity, timeMeshParams.first + 1);

    // Initialize grid (arma::mat is used for the 2D grid)
    grid = arma::zeros<arma::mat>(assetPrices.n_elem, timeSteps.n_elem);
}

// Apply boundary conditions
void Solver::applyBoundaryConditions() {
    // Terminal condition (at maturity T)
    for (size_t i = 0; i < assetPrices.n_elem; ++i) {
        if (optionType == "call")
            grid(i, timeSteps.n_elem - 1) = std::max(0.0, assetPrices(i) - strike);
        else if (optionType == "put")
            grid(i, timeSteps.n_elem - 1) = std::max(0.0, strike - assetPrices(i));
    }

    // Boundary conditions for American/European options
    for (size_t t = 0; t < timeSteps.n_elem; ++t) {
        grid(0, t) = (optionType == "call") ? 0.0 : (strike * std::exp(-interpolateRate(timeSteps(t)) * (maturity - timeSteps(t))));
        grid(assetPrices.n_elem - 1, t) = (optionType == "call") ? (assetPrices(assetPrices.n_elem - 1) - strike) : 0.0;
    }
}

// Crank-Nicholson step
void Solver::crankNicholsonStep(size_t t) {
    size_t n = assetPrices.n_elem - 2; // Exclude boundary points
    arma::vec V_prev(n);
    arma::vec V_next(n);

    // Extract the relevant values for V_prev (exclude boundaries)
    for (size_t i = 1; i <= n; ++i) {
        V_prev(i - 1) = grid(i, t);
    }

    // Define coefficients
    arma::vec alpha(n), beta(n), gamma(n);
    double dt = timeSteps(t) - timeSteps(t - 1);
    double dS = assetPrices(1) - assetPrices(0);

    for (size_t i = 0; i < n; ++i) {
        double S = assetPrices(i + 1);
        alpha(i) = 0.25 * dt * ((volatility * volatility * S * S) / (dS * dS) - interpolateRate(timeSteps(t)) * S / dS);
        beta(i) = -0.5 * dt * ((volatility * volatility * S * S) / (dS * dS) + interpolateRate(timeSteps(t)));
        gamma(i) = 0.25 * dt * ((volatility * volatility * S * S) / (dS * dS) + interpolateRate(timeSteps(t)) * S / dS);
    }

    // Create tridiagonal matrices
    arma::sp_mat ML(n, n), MR(n, n);
    ML.diag(-1) = -alpha.subvec(1, n - 1);
    ML.diag(0) = 1.0 - beta;
    ML.diag(1) = -gamma.subvec(0, n - 2);

    MR.diag(-1) = alpha.subvec(1, n - 1);
    MR.diag(0) = 1.0 + beta;
    MR.diag(1) = gamma.subvec(0, n - 2);

    // Compute the right-hand side: MR * V_prev
    V_next = MR * V_prev;

    // Solve the linear system: ML * V_prev = V_next
    arma::spsolve(V_prev, ML, V_next, "lapack");

    // Update the grid with the computed values
    for (size_t i = 1; i <= n; ++i) {
        grid(i, t - 1) = V_prev(i - 1);
    }
}

// Interpolate the risk-free rate
double Solver::interpolateRate(double time) {
    auto it = std::lower_bound(riskFreeRate.begin(), riskFreeRate.end(),
                               std::make_pair(time, 0.0),
                               [](const std::pair<double, double>& a, const std::pair<double, double>& b) { return a.first < b.first; });

    if (it == riskFreeRate.begin())
        return it->second;
    if (it == riskFreeRate.end())
        return (it - 1)->second;

    double t1 = (it - 1)->first, r1 = (it - 1)->second;
    double t2 = it->first, r2 = it->second;
    return r1 + (time - t1) * (r2 - r1) / (t2 - t1);
}

// Getter for option type
std::string Solver::getOptionType() const {
    return optionType;
}

// Getters for grid, asset prices, and time steps
const arma::mat& Solver::getGrid() const {
    return grid;
}

const arma::vec& Solver::getAssetPrices() const {
    return assetPrices;
}

const arma::vec& Solver::getTimeSteps() const {
    return timeSteps;
}

// Méthode pour calculer le prix d'un Call Européen avec Crank-Nicholson
double Solver::CrankNicolson_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    // Initialisation des paramètres
    optionType = "call";  // On spécifie que c'est un Call Européen
    strike = K;
    maturity = T;
    spotPrice = S0;
    
    // Configuration de la grille avec les paramètres donnés
    timeMeshParams = std::make_pair(dis_tmp, dis_tmp);  // Nombre de pas de temps (en général T / dt)
    spotMeshParams = std::make_pair(dis_spac, dis_spac);  // Nombre de pas pour l'actif sous-jacent (en général S_max / dS)
    
    // Reconfigurer la grille et appliquer les conditions aux bords
    setupGrid(); 
    applyBoundaryConditions();

    // Calculer le prix en utilisant la méthode Crank-Nicholson
    for (size_t t = timeSteps.n_elem - 1; t > 0; --t) {
        crankNicholsonStep(t);
    }

    // Retourner le prix de l'option Call Européenne à t=0
    return grid(spotMeshParams.first / 2, 0);  // On prend la valeur centrale pour S0
}

// Méthode pour calculer le prix d'un Put Européen avec Crank-Nicholson
double Solver::CrankNicolson_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    // Initialisation des paramètres
    optionType = "put";  // On spécifie que c'est un Put Européen
    strike = K;
    maturity = T;
    spotPrice = S0;
    
    // Configuration de la grille avec les paramètres donnés
    timeMeshParams = std::make_pair(dis_tmp, dis_tmp);  // Nombre de pas de temps (en général T / dt)
    spotMeshParams = std::make_pair(dis_spac, dis_spac);  // Nombre de pas pour l'actif sous-jacent (en général S_max / dS)
    
    // Reconfigurer la grille et appliquer les conditions aux bords
    setupGrid(); 
    applyBoundaryConditions();

    // Calculer le prix en utilisant la méthode Crank-Nicholson
    for (size_t t = timeSteps.n_elem - 1; t > 0; --t) {
        crankNicholsonStep(t);
    }

    // Retourner le prix de l'option Put Européenne à t=0
    return grid(spotMeshParams.first / 2, 0);  // On prend la valeur centrale pour S0
}

// Méthode pour calculer le prix d'un Call Américain avec Crank-Nicholson
double Solver::CrankNicolson_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac) {
    // Initialisation des paramètres
    optionType = "call";  // On spécifie que c'est un Call Américain
    strike = K;
    maturity = T;
    spotPrice = S0;
    
    // Configuration de la grille avec les paramètres donnés
    timeMeshParams = std::make_pair(dis_tmp, dis_tmp);  // Nombre de pas de temps (en général T / dt)
    spotMeshParams = std::make_pair(dis_spac, dis_spac);  // Nombre de pas pour l'actif sous-jacent (en général S_max / dS)
    
    // Reconfigurer la grille et appliquer les conditions aux bords
    setupGrid(); 
    applyBoundaryConditions();

    // Calculer le prix en utilisant la méthode Crank-Nicholson
    for (size_t t = timeSteps.n_elem - 1; t > 0; --t) {
        crankNicholsonStep(t);
        // Appliquer l'exercice anticipé pour l'option américaine
        for (size_t i = 1; i < assetPrices.n_elem - 1; ++i) {
            grid(i, t - 1) = std::max(grid(i, t - 1), assetPrices(i) - strike);  // Valeur d'exercice anticipé
        }
    }

    // Retourner le prix de l'option Call Américain à t=0
    return grid(spotMeshParams.first / 2, 0);  // On prend la valeur centrale pour S0
}
