#pragma once

#include "utils.h"

void LUDecompose(const Matrix &A, Matrix &L, Matrix &U);
Vector LUSolve(const Matrix &L, const Matrix &U, const Vector &b);

