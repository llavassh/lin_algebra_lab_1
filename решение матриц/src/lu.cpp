#include "../include/lu.h"
#include <cmath>

void LUDecompose(const Matrix &A, Matrix &L, Matrix &U) {
    int n = A.size();
    L = Matrix(n, Vector(n, 0.0));
    U = Matrix(n, Vector(n, 0.0));

    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }

            U[i][j] = A[i][j] - sum;
        }

        for (int j = i + 1; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[j][k] * U[k][i];
            }

            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }
}

Vector LUSolve(const Matrix &L, const Matrix &U, const Vector &b) {
    int n = L.size();
    Vector x(n, 0.0);
    Vector y(n, 0.0);

    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}
