#include "../include/experiments.h"
#include "../include/utils.h"
#include "../include/gauss.h"
#include "../include/lu.h"

#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>

void experiment1() {
    std::cout << "Experiment 1: Comparison of solution times\n" << std::endl;

    std::vector<int> sizes = {100, 200, 500, 1000};

    std::cout << std::left << std::setw(8) << "n"
              << std::setw(18) << "Gauss no pivot"
              << std::setw(18) << "Gauss pivot"
              << std::setw(18) << "LU total"
              << std::setw(18) << "LU decomp"
              << std::setw(18) << "LU solve" << std::endl;
    std::cout << std::string(98, '-') << std::endl;

    for (int n : sizes) {
        std::cout << "Calculation for n = " << n << "..." << std::flush;

        Matrix A = GenerateRandomMatrix(n, 42);
        Vector b = GenerateRandomRHS(n, 43);

        double time1, time2, timeLU, timeLUDecomp, timeLUSolve;

        time1 = MeasureTime([&](){ 
            GaussNoPivot(A, b); 
        });
        time2 = MeasureTime([&](){ 
            GaussWithPivot(A, b); 
        });

        Matrix L, U;
        timeLUDecomp = MeasureTime([&](){ 
            LUDecompose(A, L, U); 
        });
        timeLUSolve = MeasureTime([&](){ 
            LUSolve(L, U, b); 
        });

        timeLU = timeLUDecomp + timeLUSolve;

        std::cout << " is ready" << std::endl;
        std::cout << std::left << std::setw(8) << n
                  << std::setw(18) << std::fixed << std::setprecision(6) << time1
                  << std::setw(18) << time2
                  << std::setw(18) << timeLU
                  << std::setw(18) << timeLUDecomp
                  << std::setw(18) << timeLUSolve << std::endl;
        std::cout << std::string(98, '-') << std::endl;
    }
}

void experiment2() {
    std::cout << "Experiment 2: Multiple Right-Hand Sides\n" << std::endl;
    int n = 500;
    std::cout << "Matrix size: n =" << n << std::endl;
    Matrix A = GenerateRandomMatrix(n, 42);

    std::vector<int> k_values = {1, 10, 100};

    std::cout << std::left << std::setw(10) << "k"
              << std::setw(20) << "Gauss pivot"
              << std::setw(20) << "LU method"
              << std::setw(15) << "Speedup" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    for (int k : k_values) {
        std::cout << "Calculation for k = " << k << "..." << std::flush;

        std::vector<Vector> B;
        for (int i = 0; i < k; i++) {
            B.push_back(GenerateRandomRHS(n, 100 + i));
        }

        double timeGauss = MeasureTime([&]() {
            for (const auto &b : B) {
                GaussWithPivot(A, b);
            }
        });

        double timeLU = MeasureTime([&]() {
            Matrix L, U;
            LUDecompose(A, L, U);
            for (const auto& b : B)
                LUSolve(L, U, b);
        });

        std::cout << " is ready" << std::endl;
        std::cout << std::left << std::setw(10) << k
                  << std::setw(20) << std::fixed << std::setprecision(6) << timeGauss
                  << std::setw(20) << timeLU
                  << std::setw(15) << (timeGauss / timeLU) << std::endl;
    }

    std::cout << std::string(65, '-') << std::endl;         
}

void experiment3() {
    std::cout << "Experiment 3: Precision on the Hilbert Matrix\n";

    std::vector<int> sizes = {5, 10, 15};

    std::cout << std::left << std::setw(8) << "n"
              << std::setw(25) << "Method"
              << std::setw(20) << "Relative error"
              << std::setw(20) << "Residual" << std::endl;
    std::cout << std::string(73, '-') << std::endl;
    
    for (int n : sizes) {
        Matrix H = MatrixHilbert(n);
        Vector x_true(n, 1.0);

        Vector b(n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                b[i] += H[i][j] * x_true[j];
            }
        }

        Vector x1 = GaussNoPivot(H, b);
        std::cout << std::left << std::setw(8) << n
                  << std::setw(25) << "No pivot"
                  << std::setw(20) << std::scientific << std::setprecision(6) << RelativeError(x_true, x1)
                  << std::setw(20) << ComputeResidual(H, x1, b) << std::endl;

        Vector x2 = GaussWithPivot(H, b);
        std::cout << std::left << std::setw(8) << ""
                  << std::setw(25) << "With pivot"
                  << std::setw(20) << std::scientific << std::setprecision(6) << RelativeError(x_true, x2)
                  << std::setw(20) << ComputeResidual(H, x2, b) << std::endl;

        std::cout << std::string(73, '-') << std::endl;
    }
}