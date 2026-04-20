#pragma once

#include <vector>
#include <random>
#include <functional>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix GenerateRandomMatrix(int n, unsigned seed = 42);
Vector GenerateRandomRHS(int n, unsigned seed = 43);

Matrix MatrixHilbert(int n);

double ComputeResidual(const Matrix &A, const Vector &x, const Vector &b);
double RelativeError(const Vector &x_true, const Vector &x_approx);

double MeasureTime(const std::function<void()> &func);