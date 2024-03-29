// Drupal : Check if the operating system is not Windows
#ifdef _WIN32
#pragma once

#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>

using namespace std;

// Variables for blas
static char chy = 'T';
static char chn = 'N';
static double one = 1.0;
static double zero = 0.0;
static int oneI = 1;

static const char* jobzV = (char*)"V";
static const char* uploU = (char*)"U";
// static const double eigenTol = -1;
static int info;
static const int workSize = 102;

using namespace std;

//JD
inline void FEMTECH_dgemv_(char trans, int m, int n, double alpha,
	const double* A, int lda,
	const double* x, int incx,
	double beta, double* y, int incy) {

	// Perform matrix-vector multiplication
	if (trans == 'N')
	{
		for (int i = 0; i < m; ++i)
		{
			double sum = 0.0;
			for (int j = 0; j < n; ++j)
			{
				sum += A[i * lda + j] * x[j * incx];
			}
			y[i * incy] = alpha * sum + beta * y[i * incy];
		}
	}
	else if (trans == 'T') {
		// Similar logic for transpose case
		for (int j = 0; j < n; ++j)
		{
			double temp = alpha * x[j * incx];
			for (int i = 0; i < m; ++i)
			{
				y[i * incy] += temp * A[i * lda + j];
			}
		}
	}
	else
	{
		/*std::cerr << "Invalid value for trans. Should be 'N' or 'T'." << std::endl;*/
		FILE_LOG_MASTER(ERROR, "Invalid value for trans. Should be 'N' or 'T'.");
		return;
	}

	std::cout << "FEMTECH_dgemv_ result:" << std::endl;
	for (int i = 0; i < m; ++i)
	{
		std::cout << y[i] << " ";
	}


}

extern "C"
{
	inline void FEMTECH_dgesv_(const int* const n, const int* const nrhs,
		double* const a, const int* const lda, int* const ipiv,
		double* const b, const int* const ldb, int* const info)
	{
		// LU decomposition with partial pivoting
		*info = 0;
		for (int k = 0; k < *n - 1; ++k)
		{
			// Find the pivot
			int pivot = k;
			for (int i = k + 1; i < *n; ++i)
			{
				if (std::abs(a[i * *lda + k]) > std::abs(a[pivot * *lda + k]))
				{
					pivot = i;
				}
			}

			// Swap rows in A and elements in ipiv
			if (pivot != k)
			{
				std::swap(ipiv[k], ipiv[pivot]);
				for (int j = 0; j < *n; ++j)
				{
					std::swap(a[k * *lda + j], a[pivot * *lda + j]);
				}
			}

			// Factorize
			if (a[k * *lda + k] == 0.0)
			{
				*info = k + 1;
				return;
			}

			for (int i = k + 1; i < *n; ++i)
			{
				a[i * *lda + k] /= a[k * *lda + k];
				for (int j = k + 1; j < *n; ++j)
				{
					a[i * *lda + j] -= a[i * *lda + k] * a[k * *lda + j];
				}
			}
		}

		// Solve the system
		for (int i = 0; i < *nrhs; ++i)
		{
			// Apply pivots to the right-hand side
			for (int k = 0; k < *n - 1; ++k)
			{
				std::swap(b[k * *ldb + i], b[ipiv[k] * *ldb + i]);
			}

			// Solve L * y = b
			for (int k = 0; k < *n - 1; ++k)
			{
				for (int j = k + 1; j < *n; ++j)
				{
					b[j * *ldb + i] -= a[j * *lda + k] * b[k * *ldb + i];
				}
			}

			// Solve U * x = y
			for (int k = *n - 1; k >= 0; --k)
			{
				b[k * *ldb + i] /= a[k * *lda + k];
				for (int j = 0; j < k; ++j)
				{
					b[j * *ldb + i] -= a[j * *lda + k] * b[k * *ldb + i];
				}
			}
		}
	}


	inline void FEMTECH_dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info)
	{
		// LU decomposition with partial pivoting
		*info = 0;

		for (int k = 0; k < std::min(*m, *n); ++k)
		{
			// Find the pivot
			int pivot = k;
			for (int i = k + 1; i < *m; ++i)
			{
				if (std::abs(a[i * *lda + k]) > std::abs(a[pivot * *lda + k]))
				{
					pivot = i;
				}
			}

			// Swap rows in A and elements in ipiv
			if (pivot != k)
			{
				std::swap(ipiv[k], ipiv[pivot]);
				for (int j = 0; j < *n; ++j)
				{
					std::swap(a[k * *lda + j], a[pivot * *lda + j]);
				}
			}

			// Factorize
			if (a[k * *lda + k] == 0.0)
			{
				*info = k + 1;
				return;
			}

			for (int i = k + 1; i < *m; ++i)
			{
				a[i * *lda + k] /= a[k * *lda + k];
				for (int j = k + 1; j < *n; ++j)
				{
					a[i * *lda + j] -= a[i * *lda + k] * a[k * *lda + j];
				}
			}
		}
	}

	inline void FEMTECH_dgetrs_(char trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info)
	{
		// Solve the system using the LU-factored matrix
		*info = 0;

		// Apply pivots to the right-hand side
		for (int k = 0; k < *n - 1; ++k)
		{
			std::swap(b[k], b[ipiv[k]]);
		}

		// Solve L * y = b
		for (int k = 0; k < *n - 1; ++k)
		{
			for (int j = k + 1; j < *n; ++j)
			{
				b[j] -= a[j * *lda + k] * b[k];
			}
		}

		// Solve U * x = y
		for (int k = *n - 1; k >= 0; --k)
		{
			b[k] /= a[k * *lda + k];
			for (int j = 0; j < k; ++j)
			{
				b[j] -= a[j * *lda + k] * b[k];
			}
		}
	}
	inline void FEMTECH_dger_(const int* const M, const int* const N, const double* const alpha, const double* const X,
		const int* const incX, const double* const Y, const int* const incY, double* const A,
		const int* const lda)
	{
		// Perform the rank-1 update
		for (int j = 0; j < *N; ++j)
		{
			for (int i = 0; i < *M; ++i)
			{
				A[i * *lda + j] += (*alpha) * X[i * *incX] * Y[j * *incY];
			}
		}
	}


	// Function to perform Jacobi rotation
	inline void jacobi_rotation(std::vector<std::vector<double>>& A, int p, int q, int n, std::vector<std::vector<double>>& V) {
		double theta = atan2(A[q][q] - A[p][p], 2 * A[p][q]);
		double c = cos(theta);
		double s = sin(theta);
		double t;

		// Update elements of A
		A[p][p] += t = c * c * (A[p][p] - A[q][q]) + s * s * (A[q][q] - A[p][p]) + 2 * c * s * A[p][q];
		A[q][q] -= t;
		A[p][q] = A[q][p] = 0;

		// Update elements of V (eigenvectors)
		for (int i = 0; i < n; ++i) {
			if (i != p && i != q) {
				A[i][p] = A[p][i] = c * A[i][p] - s * A[i][q];
				A[i][q] = A[q][i] = c * A[i][q] + s * A[i][p];
			}
			t = c * V[i][p] - s * V[i][q];
			V[i][q] = s * V[i][p] + c * V[i][q];
			V[i][p] = t;
		}

		// Print Jacobi rotation step
		/*std::cout << "Jacobi Rotation: (p, q) = (" << p << ", " << q << "), theta = " << theta << std::endl;*/
	}

	// Function to perform Jacobi eigenvalue algorithm
	inline void jacobi_eigenvalue_algorithm(std::vector<std::vector<double>>& A, int n, std::vector<double>& eigenvalues, std::vector<std::vector<double>>& eigenvectors) {
		eigenvectors = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
		for (int i = 0; i < n; ++i) eigenvectors[i][i] = 1;
		int max_iterations = 1000;
		double tolerance = 1e-10;
		int iterations = 0;
		double off_diagonal_norm = 1.0;

		while (off_diagonal_norm > tolerance && iterations < max_iterations) {
			off_diagonal_norm = 0.0;
			for (int i = 0; i < n - 1; ++i) {
				for (int j = i + 1; j < n; ++j) {
					off_diagonal_norm += A[i][j] * A[i][j];
				}
			}
			off_diagonal_norm = std::sqrt(off_diagonal_norm);

			if (off_diagonal_norm <= tolerance) break;

			for (int i = 0; i < n - 1; ++i) {
				for (int j = i + 1; j < n; ++j) {
					if (std::abs(A[i][j]) > tolerance * off_diagonal_norm) {
						jacobi_rotation(A, i, j, n, eigenvectors);
					}
				}
			}

			++iterations;
		}

		for (int i = 0; i < n; ++i) eigenvalues[i] = A[i][i];
	}

	// FEMTECH_dsyev_ declaration
	inline void FEMTECH_dsyev_(const char* const jobz, const char* const uplo,
		const int* const n, double* const a, const int* const lda,
		double* const w, double* work, const int* const lwork, int* const info) {
		// Check for valid parameters
		if (*n <= 0 || *lda < *n || *lwork < *n) {
			*info = -3;
			return;
		}

		// Make a copy of the input matrix
		std::vector<std::vector<double>> A(*n, std::vector<double>(*n, 0));
		for (int i = 0; i < *n; ++i) {
			for (int j = 0; j < *n; ++j) {
				A[i][j] = a[i * *lda + j];
			}
		}

		// Compute eigenvalues and eigenvectors using Jacobi eigenvalue algorithm
		std::vector<double> eigenvalues(*n);
		std::vector<std::vector<double>> eigenvectors;
		jacobi_eigenvalue_algorithm(A, *n, eigenvalues, eigenvectors);

		//// Copy eigenvalues to output array
		//for (int i = 0; i < *n; ++i) {
		//	w[i] = eigenvalues[i];
		//}

		std::copy(eigenvalues.begin(), eigenvalues.end(), w);

		// Set info to indicate success
		*info = 0;
	}
}

// Perform matrix-matrix multiplication
inline void FEMTECH_dgemm_(char transa, char transb, int m, int n, int k, double alpha,
	const double* A, int lda,
	const double* B, int ldb, double beta,
	double* C, int ldc) {

	if (transa == 'N' && transb == 'N') {
		// No transpose for both A and B
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				double sum = 0.0;
				for (int l = 0; l < k; ++l) {
					sum += A[i * lda + l] * B[l * ldb + j];
				}
				C[i * ldc + j] = alpha * sum + beta * C[i * ldc + j];
			}
		}
	}
	else if (transa == 'T' && transb == 'N') {
		// Transpose A, No transpose B
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				double sum = 0.0;
				for (int l = 0; l < k; ++l) {
					sum += A[l * lda + i] * B[l * ldb + j];
				}
				C[i * ldc + j] = alpha * sum + beta * C[i * ldc + j];
			}
		}
	}
	else if (transa == 'N' && transb == 'T') {
		// No transpose A, Transpose B
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				double sum = 0.0;
				for (int l = 0; l < k; ++l) {
					sum += A[i * lda + l] * B[j * ldb + l];
				}
				C[i * ldc + j] = alpha * sum + beta * C[i * ldc + j];
			}
		}
	}
	else if (transa == 'T' && transb == 'T') {
		// Transpose A, Transpose B
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				double sum = 0.0;
				for (int l = 0; l < k; ++l) {
					sum += A[l * lda + i] * B[j * ldb + l];
				}
				C[i * ldc + j] = alpha * sum + beta * C[i * ldc + j];
			}
		}
	}
	else {
		std::cerr << "Invalid transpose options." << std::endl;
	}
}

#endif // MATRIX_OPERATIONS_H
#endif // _WIN32