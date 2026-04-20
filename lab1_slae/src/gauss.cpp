#include "../include/gauss.h"
#include <algorithm>
#include <cmath>

Vector GaussNoPivot(const Matrix &A, const Vector &b) {
    int n = A.size();
    Matrix augmented(n, Vector(n + 1));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n] = b[i];
    }

    for (int k = 0; k < n - 1; k++) {
        if (std::abs(augmented[k][k]) < 1e-12) {
            return Vector(n, 0.0);
        }

        for (int i = k + 1; i < n; i++) {
            double factor = augmented[i][k] / augmented[k][k];
            for (int j = k; j <= n; j++) {
            augmented[i][j] -= factor * augmented[k][j];
            }
        }
    }

    Vector x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += augmented[i][j] * x[j];
        }

        x[i] = (augmented[i][n] - sum) / augmented[i][i];
    }

    return x;
}

Vector GaussWithPivot(const Matrix &A, const Vector &b) {
    int n = A.size();
    Matrix augmented(n, Vector(n + 1));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n] = b[i];
    }

    for (int k = 0; k < n-1; k++) {
        int MaxRow = k;
        double MaxValue = std::abs(augmented[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (std::abs(augmented[i][k]) > MaxValue) {
                MaxValue = std::abs(augmented[i][k]);
                MaxRow = i;
            }
        }

        if (MaxRow != k) {
            std::swap(augmented[k], augmented[MaxRow]);
        }

        if (std::abs(augmented[k][k]) < 1e-12) continue;

        for (int i = k + 1; i < n; i++) {
            double factor = augmented[i][k] / augmented[k][k];
            for (int j = k; j <= n; j++) {
                augmented[i][j] -= factor * augmented[k][j];
            }
        }
    }

    Vector x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += augmented[i][j] * x[j];
        }

        x[i] = (augmented[i][n] - sum) / augmented[i][i];
    }

    return x;
}