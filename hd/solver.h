#pragma once

#include <string>
#include <armadillo> // Inclure Armadillo pour les calculs matriciels

class Solver {
public:
    // Constructeur
    Solver(const std::string& optionType, const std::string& exerciseType,
           double maturity, double strike, double computationDate,
           const std::pair<double, double>& timeMeshParams,
           const std::pair<double, double>& spotMeshParams,
           double spotPrice, const std::vector<std::pair<double, double>>& riskFreeRate,
           double volatility);

    // Méthode solve pour calculer les prix des options
    void solve();

    // Getters pour les résultats
    const arma::mat& getGrid() const;
    const arma::vec& getAssetPrices() const;
    const arma::vec& getTimeSteps() const;
    std::string getOptionType() const;  // Ajouter cette ligne pour récupérer le type d'option
    // Méthodes pour calculer les prix des options avec Crank-Nicholson
     double CrankNicolson_Call_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
     double CrankNicolson_Put_European(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
     double  CrankNicolson_Call_American(double S0, double K, double T, double r, double sigma, int dis_tmp, int dis_spac);
    

private:
    // Paramètres d'entrée
    std::string optionType;                   // "call" ou "put"
    std::string exerciseType;                 // "european" ou "american"
    double maturity;                          // Maturité de l'option (T)
    double strike;                            // Prix d'exercice (K)
    double computationDate;                   // Date de calcul actuelle (T0)
    double spotPrice;                         // Prix actuel de l'actif sous-jacent (S0)
    std::pair<double, double> timeMeshParams; // Paramètres de maillage temporel: (nombre de pas, taille du pas temporel)
    std::pair<double, double> spotMeshParams; // Paramètres de maillage des prix : (nombre de pas, taille du pas de prix)
    std::vector<std::pair<double, double>> riskFreeRate; // Taux sans risque sous forme de fonction linéaire par morceaux
    double volatility;                        // Volatilité (sigma)

    // Données de la grille et de calcul
    arma::mat grid;    // Grille de valeurs d'option
    arma::vec assetPrices;  // Prix discrets des actifs sous-jacents
    arma::vec timeSteps;    // Pas de temps discrets

    // Méthodes internes
    void setupGrid();                         // Initialisation de la grille et des paramètres
    void applyBoundaryConditions();           // Appliquer les conditions aux frontières de la grille
    void crankNicholsonStep(size_t timeStep); // Effectuer une étape Crank-Nicholson
    double interpolateRate(double time);      // Interpoler le taux sans risque à un instant donné
};
