
#include <cxxtest/TestSuite.h>
#include "../hd/solver.h"

class SolverTests : public CxxTest::TestSuite {
public:

    void testEuropeanCall() {
        // Arrange
        Solver solver("call", "european", 1.0, 100.0, 0.0, {100, 0.01}, {100, 1.0},
                      100.0, {{0.0, 0.05}, {1.0, 0.05}}, 0.2);

        // Act
        solver.solve();
        const auto& grid = solver.getGrid();
        const auto& assetPrices = solver.getAssetPrices();

        // Assert
        double spotPrice = 100.0;
        auto it = std::find(assetPrices.begin(), assetPrices.end(), spotPrice);
        TS_ASSERT(it != assetPrices.end());

        size_t idx = std::distance(assetPrices.begin(), it);
        double priceAtSpot = grid[idx][0];

        // Known analytical price using Black-Scholes formula
        double analyticalPrice = 10.45058; // Replace with the actual computed value
        TS_ASSERT_DELTA(priceAtSpot, analyticalPrice, 1e-2); // Allow 1% tolerance
    }

    void testEuropeanPut() {
        // Arrange
        Solver solver("put", "european", 1.0, 100.0, 0.0, {100, 0.01}, {100, 1.0},
                      100.0, {{0.0, 0.05}, {1.0, 0.05}}, 0.2);

        // Act
        solver.solve();
        const auto& grid = solver.getGrid();
        const auto& assetPrices = solver.getAssetPrices();

        // Assert
        double spotPrice = 100.0;
        auto it = std::find(assetPrices.begin(), assetPrices.end(), spotPrice);
        TS_ASSERT(it != assetPrices.end());

        size_t idx = std::distance(assetPrices.begin(), it);
        double priceAtSpot = grid[idx][0];

        // Known analytical price using Black-Scholes formula
        double analyticalPrice = 5.57352; // Replace with the actual computed value
        TS_ASSERT_DELTA(priceAtSpot, analyticalPrice, 1e-2); // Allow 1% tolerance
    }

    void testBoundaryConditions() {
        // Arrange
        Solver solver("call", "european", 1.0, 100.0, 0.0, {100, 0.01}, {100, 1.0},
                      100.0, {{0.0, 0.05}, {1.0, 0.05}}, 0.2);

        // Act
        solver.solve();
        const auto& grid = solver.getGrid();
        const auto& assetPrices = solver.getAssetPrices();

        // Assert: Check boundary values
        TS_ASSERT_EQUALS(grid.front().front(), 0.0); // Price at S_min is 0 for call
        TS_ASSERT_DELTA(grid.back().front(), assetPrices.back() - 100.0, 1e-2); // Price at S_max
    }

    void testAmericanVsEuropean() {
        // Arrange
        Solver europeanSolver("put", "european", 1.0, 100.0, 0.0, {100, 0.01}, {100, 1.0},
                              100.0, {{0.0, 0.05}, {1.0, 0.05}}, 0.2);
        Solver americanSolver("put", "american", 1.0, 100.0, 0.0, {100, 0.01}, {100, 1.0},
                              100.0, {{0.0, 0.05}, {1.0, 0.05}}, 0.2);

        // Act
        europeanSolver.solve();
        americanSolver.solve();

        const auto& europeanGrid = europeanSolver.getGrid();
        const auto& americanGrid = americanSolver.getGrid();
        const auto& assetPrices = europeanSolver.getAssetPrices();

        // Assert: American price should not be less than European price
        for (size_t i = 0; i < assetPrices.size(); ++i) {
            TS_ASSERT(americanGrid[i][0] >= europeanGrid[i][0]);
        }
    }
};
