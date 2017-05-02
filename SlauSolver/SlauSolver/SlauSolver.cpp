#include "stdafx.h"
#include "SlauSolver.h"

CSlauSolver::CSlauSolver()
{
	return;
}

double CSlauSolver::Dot(double *a, double *b, int n)
{
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += a[i] * b[i];
	}

	return res;
}

void CSlauSolver::Sum(double *a, double *b, int n, double alfa)
{
	for (int i = 0; i < n; i++)
	{
		a[i] += b[i] * alfa;
	}
}

void CSlauSolver::Sum(double *a, double *b, double * res, int n, double alfa)
{
	for (int i = 0; i < n; i++)
	{
		res[i] = a[i] + b[i] * alfa;
	}
}


void CSlauSolver::Diff(double *a, double *b, int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] -= b[i];
	}
}

void CSlauSolver::Mult(CRSMatrix & A, double * b, double * res)
{
	int elemCountInRow = 0;

	for (int i = 0; i < A.n; i++)
	{
		elemCountInRow = A.rowPtr[i + 1] - A.rowPtr[i];
		res[i] = 0;

		for (int j = 0; j < elemCountInRow; j++)
		{
			res[i] += A.val[A.rowPtr[i] + j] * b[A.colIndex[A.rowPtr[i] + j]];
		}
	}
}

void CSlauSolver::SolveR(CRSMatrix & A, double * z, double * r, double alfa)
{
	int elemCountInRow = 0;
	double res;
	for (int i = 0; i < A.n; i++)
	{
		elemCountInRow = A.rowPtr[i + 1] - A.rowPtr[i];
		res = 0;

		for (int j = 0; j < elemCountInRow; j++)
		{
			res += A.val[A.rowPtr[i] + j] * z[A.colIndex[A.rowPtr[i] + j]] * alfa;
		}

		r[i] += res;
	}
}

void  CSlauSolver::GenerateSolution(double * x, int n)
{
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}
}

void CSlauSolver::Copy(double * a, double * copyA, int n)
{
	for (int i = 0; i < n; i++)
	{
		copyA[i] = a[i];
	}
}

bool CSlauSolver::IsEnd(double * x0, double * x, int n, double eps)
{
	for (int i = 0; i < n; i++)
	{
		if (abs(x0[i] - x[i]) > eps)
		{
			return false;
		}
	}

	return true;
}

void CSlauSolver::SLE_Solver_CRS_BICG(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{
	int n = A.n;

	double * r = new double[n]();
	double * p = new double[n]();
	double * z = new double[n]();
	double * s = new double[n]();
	double * x0 = new double[n]();
	double * temp = new double[n]();

	double * predR = new double[n]();
	double * predP = new double[n]();

	double alfa, betta;


	GenerateSolution(x, n);
	SolveR(A, x, r, -1);
	Copy(r, p, n);
	Copy(r, z, n);
	Copy(r, s, n);
	Copy(x, x0, n);

	while (!IsEnd(x0, x, n, eps))
	{
		Copy(x, x0, n);
		Copy(r, predR, n);
		Copy(p, predP, n);

		Mult(A, z, temp);
		alfa = Dot(p, r, n) / Dot(s, temp, n);

		Sum(x, z, n, alfa); // пересчет x

		SolveR(A, z, r, -alfa); // пересчет r

		SolveR(A, s, p, -alfa); // пересчет p

		betta = Dot(p, r, n) / Dot(predP, predR, n);

		Sum(r, z, z, betta); // пересчет z

		Sum(p, s, s, betta); // пересчет s
	}

	delete[] r, p, z, s, x, temp;
}




