#include "stdafx.h"
#include "SlauSolver.h"

CSlauSolver::CSlauSolver()
{
	return;
}

double CSlauSolver::Dot(double *a, double *b, int n)
{
	double res = 0;
//#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		res += a[i] * b[i];
	}

	return res;
}

void CSlauSolver::Sum(double *a, double *b, int n, double alfa)
{
//#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		a[i] += b[i] * alfa;
	}
}

void CSlauSolver::Sum(double *a, double *b, double * res, int n, double alfa)
{
//#pragma omp parallel for
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

//#pragma omp parallel for
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

void CSlauSolver::SolveR(CRSMatrix & A, double * z, double * b, double * r, double alfa)
{
	int elemCountInRow = 0;
	double res;

//#pragma omp parallel for
	for (int i = 0; i < A.n; i++)
	{
		elemCountInRow = A.rowPtr[i + 1] - A.rowPtr[i];
		res = 0;

		for (int j = 0; j < elemCountInRow; j++)
		{
			res += A.val[A.rowPtr[i] + j] * z[A.colIndex[A.rowPtr[i] + j]] * alfa;
		}

		r[i] = b[i] + res;
	}
}

void CSlauSolver::SolveRT(CRSMatrix & A, double * z, double * b, double * r, double alfa)
{
	int elemCountInCol = 0;
	double res;

//#pragma omp parallel for
	for (int i = 0; i < A.m; i++)
	{
		res = 0;
		for (int j = 0; j < A.colIndex.size(); j++)
		{
			if (A.colIndex[j] == i)
			{
				res += alfa * A.val[j] * z[GetRowIndex(A, j)];
			}
		}

		r[i] = b[i] + res;
	}
}

int CSlauSolver::GetRowIndex(CRSMatrix & A, int index)
{
	int res;
	for (int i = 0; i < A.rowPtr.size() - 1; i++)
	{
		if (index >= A.rowPtr[i] && index < A.rowPtr[i + 1])
		{
			res = i;
			break;
		}
	}

	return res;
}

void  CSlauSolver::GenerateSolution(double * x, int n)
{
	for (int i = 0; i < n; i++)
	{
		x[i] = 1;
	}
}

void CSlauSolver::Copy(double * a, double * copyA, int n)
{
	//#pragma omp parallel for
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
	SolveR(A, x, b, r, -1);
	Copy(r, p, n);
	Copy(r, z, n);
	Copy(r, s, n);

	do
	{
		count++;
		Copy(x, x0, n);
		Copy(r, predR, n);
		Copy(p, predP, n);

		Mult(A, z, temp);
		alfa = Dot(p, r, n) / Dot(s, temp, n);

		Sum(x, z, n, alfa); // пересчет x

		SolveR(A, z, r, r, -alfa); // пересчет r

		SolveRT(A, s, p, p, -alfa); // пересчет p

		betta = Dot(p, r, n) / Dot(predP, predR, n);

		Sum(r, z, z, n, betta); // пересчет z

		Sum(p, s, s, n, betta); // пересчет s

		if (count >= max_iter)
			break;

	} while (!IsEnd(x0, x, n, eps));

	delete[] r, p, z, s, x, x0, temp, predR, predP;
}



