#include "../hd/solver.h"
#include "../hd/greeks.h"
#include <iostream>
#include <vector>
#include <utility>  // Pour pair

int main() {
    // Paramètres de test
    std::string optionType = "call";           // Type d'option : "call" ou "put"
    std::string exerciseType = "European";     // Type d'exercice : "European" ou "American"
    double maturity = 1.0;                     // Maturité de l'option (en années)
    double strike = 100.0;                     // Strike price
    double computationDate = 0.0;              // Date de computation (généralement 0)
    
    // Paramètres du maillage en temps et en prix de l'actif sous-jacent
    std::pair<double, double> timeMeshParams(100, 1.0);  // 100 pas de temps, maturité 1.0
    std::pair<double, double> spotMeshParams(100, 1.0);  // 100 pas de prix, prix initial 1.0

    double spotPrice = 100.0;  // Prix initial de l'actif sous-jacent

    // Taux d'intérêt sans risque (exemple : de 0 à 1 an avec taux constant)
    std::vector<std::pair<double, double>> riskFreeRate = { {0.0, 0.05}, {1.0, 0.05} };

    double volatility = 0.2;  // Volatilité

    // Créer l'objet Solver
    Solver solver(optionType, exerciseType, maturity, strike, computationDate,
                  timeMeshParams, spotMeshParams, spotPrice, riskFreeRate, volatility);
    
    // Résoudre le problème
    solver.solve();

    // Afficher les résultats de la grille (grille de valeurs de l'option)
    const auto& grid = solver.getGrid();
    const auto& assetPrices = solver.getAssetPrices();
    const auto& timeSteps = solver.getTimeSteps();

    std::cout << "Grid values at maturity (last time step):\n";
    for (size_t i = 0; i < assetPrices.size(); ++i) {
        std::cout << "Asset price: " << assetPrices[i]
                  << ", Option value: " << grid(i, grid.n_cols - 1) << std::endl;  // Accéder à la dernière colonne de chaque ligne
    }

    std::cout << "\nTime steps:\n";
    for (size_t t = 0; t < timeSteps.size(); ++t) {
        std::cout << "Time step " << t << ": " << timeSteps[t] << " years\n";
    }

    // Créer l'objet Greeks
    Greeks greeks(solver);  // Crée une instance de la classe Greeks

    // Calcul des Greeks
    if (optionType == "call" && exerciseType == "European") {
        // Calcul des Greeks pour une option Call Européenne
        double delta = greeks.Delta_Call_European(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);
        double gamma = greeks.Gamma_Call_European(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);
        double theta = greeks.Theta_Call_European(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);
        double vega = greeks.Vega_Call_European(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);
        double rho = greeks.Rho_Call_European(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);

        // Afficher les Greeks calculés pour l'option Call Européenne
        std::cout << "\nGreeks for European Call option:\n";
        std::cout << "Delta: " << delta << "\n";
        std::cout << "Gamma: " << gamma << "\n";
        std::cout << "Theta: " << theta << "\n";
        std::cout << "Vega: " << vega << "\n";
        std::cout << "Rho: " << rho << "\n";
    } else if (optionType == "call" && exerciseType == "American") {
        // Calcul des Greeks pour une option Call Américaine
        double delta = greeks.Delta_Call_American(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);
        double gamma = greeks.Gamma_Call_American(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);
        double theta = greeks.Theta_Call_American(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);
        double vega = greeks.Vega_Call_American(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);
        double rho = greeks.Rho_Call_American(spotPrice, strike, maturity, riskFreeRate[0].second, volatility, timeMeshParams.first, spotMeshParams.first);

        // Afficher les Greeks calculés pour l'option Call Américaine
        std::cout << "\nGreeks for American Call option:\n";
        std::cout << "Delta: " << delta << "\n";
        std::cout << "Gamma: " << gamma << "\n";
        std::cout << "Theta: " << theta << "\n";
        std::cout << "Vega: " << vega << "\n";
        std::cout << "Rho: " << rho << "\n";
    }

    return 0;
}
