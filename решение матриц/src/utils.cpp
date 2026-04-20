#include "../include/utils.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <random>
#include <functional>

Matrix GenerateRandomMatrix(int n, unsigned seed) {
    if (n == 0) {
        return Matrix{};
    }

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    Matrix A(n, Vector(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = dist(rng);
        }
    }

    return A;
}

Vector GenerateRandomRHS(int n, unsigned seed) {
    if (n == 0) {
        return Vector{};
    }

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    Vector b(n);
    for (int i = 0; i < n; i++) {
        b[i] = dist(rng);
    }

    return b;
}

Matrix MatrixHilbert(int n) {
    Matrix H(n, Vector(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }

    return H;
}

double ComputeResidual(const Matrix &A, const Vector &x, const Vector &b) {
    int n = A.size();
    Vector Ax(n, 0.0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Ax[i] += A[i][j] * x[j];
        }
    }

    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = Ax[i] - b[i];
        norm += diff * diff;
    }

    return std::sqrt(norm);
}

double RelativeError(const Vector &x_true, const Vector &x_approx) {
    int n = x_true.size();
    double norm_true = 0.0, norm_diff = 0.0;

    for (int i = 0; i < n; i++) {
        norm_true += x_true[i] * x_true[i];
        double diff = x_true[i] - x_approx[i];
        norm_diff += diff * diff; 
    }

    return std::sqrt(norm_diff/norm_true);
}

double MeasureTime(const std::function<void()>& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count();
}