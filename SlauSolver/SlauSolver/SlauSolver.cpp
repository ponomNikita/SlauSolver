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

bool CSlauSolver::IsEnd(double * x, int n, double eps)
{
	double norma = 0;
	for (int i = 0; i < n; i++)
	{
		norma += x[i] * x[i];
	}

	norma = sqrt(norma);
	return norma < eps;
}

void CSlauSolver::SLE_Solver_CRS_BICG(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{
	int n = A.n;

	double * r = new double[n]();
	double * p = new double[n]();
	double * r_sop = new double[n]();
	double * p_sop = new double[n]();
	double * temp = new double[n]();
	double * predR_sop = new double[n]();
	double * predR = new double[n]();

	double alfa, betta;


	GenerateSolution(x, n); // начальное приближение
	SolveR(A, x, b, r, -1); // начальное r
	Copy(r, p, n);
	Copy(r, r_sop, n);
	Copy(r, p_sop, n);

	do
	{
		count++;

		Copy(r, predR, n);
		Copy(r_sop, predR_sop, n);

		Mult(A, p, temp);
		alfa = Dot(r, r_sop, n) / Dot(temp, p_sop, n);

		Sum(x, p, x, n, alfa);

		SolveR(A, p, r, r, -alfa);
		SolveRT(A, p_sop, r_sop, r_sop, -alfa);

		betta = Dot(r, r_sop, n) / Dot(predR, predR_sop, n);

		if (betta == 0 || IsEnd(r, n, eps))
			break;

		if (count >= max_iter)
			break;

		Sum(r, p, p, n, betta); // пересчт p
		Sum(r_sop, p_sop, p_sop, n, betta); // пересчт p_sop

	} while (true);

	delete[] r, p, temp, predR, r_sop, p_sop, predR_sop;
}



