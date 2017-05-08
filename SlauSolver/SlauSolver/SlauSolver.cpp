#include "stdafx.h"
#include "SlauSolver.h"
#include <omp.h>

CSlauSolver::CSlauSolver()
{
	return;
}

double CSlauSolver::Dot(double *a, double *b, int n)
{
	double res = 0;
#pragma omp parallel for reduction(+:res)
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
#pragma omp parallel for
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
#pragma omp parallel for schedule(runtime) 
	for (int i = 0; i < A.n; i++)
	{
		int elemCountInRow = A.rowPtr[i + 1] - A.rowPtr[i];
		res[i] = 0;

		for (int j = 0; j < elemCountInRow; j++)
		{
			res[i] += A.val[A.rowPtr[i] + j] * b[A.colIndex[A.rowPtr[i] + j]];
		}
	}
}

void CSlauSolver::MultTransp(CRSMatrix & A, double * b, double * res)
{
#pragma omp parallel for
	for (int i = 0; i < A.m; i++)
	{
		res[i] = 0;
		for (int j = 0; j < A.colIndex.size(); j++)
		{
			if (A.colIndex[j] == i)
			{
				res[i] += A.val[j] * b[GetRowIndex(A, j)];
			}
		}
	}
}

void CSlauSolver::SolveR(CRSMatrix & A, double * z, double * b, double * r, double alfa)
{
#pragma omp parallel for
	for (int i = 0; i < A.n; i++)
	{
		int elemCountInRow = A.rowPtr[i + 1] - A.rowPtr[i];
		double res = 0;

		for (int j = 0; j < elemCountInRow; j++)
		{
			res += A.val[A.rowPtr[i] + j] * z[A.colIndex[A.rowPtr[i] + j]] * alfa;
		}

		r[i] = b[i] + res;
	}
}

void CSlauSolver::SolveRWithResolveX(CRSMatrix & A, double * z, double * b, double * r, double * x, double * p, double alfa)
{
#pragma omp parallel for
	for (int i = 0; i < A.n; i++)
	{
		int elemCountInRow = A.rowPtr[i + 1] - A.rowPtr[i];
		double res = 0;

		for (int j = 0; j < elemCountInRow; j++)
		{
			res += A.val[A.rowPtr[i] + j] * z[A.colIndex[A.rowPtr[i] + j]] * alfa;
		}

		r[i] = b[i] + res;
		x[i] = x[i] - p[i] * alfa;
	}
}

void CSlauSolver::SolveRT(CRSMatrix & A, double * z, double * b, double * r, double alfa)
{
#pragma omp parallel for
	for (int i = 0; i < A.m; i++)
	{
		double res = 0;
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
	#pragma omp parallel for
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

double CSlauSolver::GetAlfaAndCopyPredArrays(CRSMatrix & A, double * r, double * r_sop, double * p, double * p_sop, double * temp, double * predR, double * predR_sop, int n)
{
	double alfa1 = 0;
	double alfa2 = 0;
#pragma omp parallel for reduction(+:alfa1, alfa2) 
	for (int i = 0; i < A.n; i++)
	{
		predR[i] = r[i];
		predR_sop[i] = r_sop[i];

		int elemCountInRow = A.rowPtr[i + 1] - A.rowPtr[i];
		temp[i] = 0;

		for (int j = 0; j < elemCountInRow; j++)
		{
			temp[i] += A.val[A.rowPtr[i] + j] * p[A.colIndex[A.rowPtr[i] + j]];
		}

		alfa1 += (r[i] * r_sop[i]);
		alfa2 += (temp[i] * p_sop[i]);
	}

	return alfa1 / alfa2;
}

double CSlauSolver::GetBetta(double * r, double * r_sop, double * predR, double * predR_sop, int n)
{
	double betta1 = 0;
	double betta2 = 0;

#pragma omp parallel for reduction(+:betta1, betta2)
	for (int i = 0; i < n; i++)
	{
		betta1 += r[i] * r_sop[i];
		betta2 += predR[i] * predR_sop[i];
	}

	return betta1 / betta2;
}

void CSlauSolver::ResolvePandPSop(double * p, double * p_sop, double * r, double * r_sop, int n, double betta)
{
#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		p[i] = r[i] + p[i] * betta;
		p_sop[i] = r_sop[i] + p_sop[i] * betta;
	}
}

void CSlauSolver::SLE_Solver_CRS_BICG(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{
	int n = A.n;

	double * r = new double[n]();
	double * p = new double[n]();
	double * r_sop = new double[n]();
	double * p_sop = new double[n]();
	double * a_p = new double[n]();
	double * aT_p = new double[n]();
	double * predR_sop = new double[n]();
	double * predR = new double[n]();

	double alfa, betta;


	GenerateSolution(x, n); // начальное приближение
	SolveR(A, x, b, r, -1); // начальное r

	Copy(r, p, n);
	Copy(r, r_sop, n);
	Copy(r, p_sop, n);
	Copy(r, predR_sop, n);
	Copy(r, predR, n);

	count = 0;

	do
	{
		count++;

		double alfa1 = 0;
		double alfa2 = 0;

		// Вычисление alfa
#pragma omp parallel for reduction(+:alfa1, alfa2)
		for (int i = 0; i < A.n; i++)
		{
			int elemCountInRow = A.rowPtr[i + 1] - A.rowPtr[i];
			a_p[i] = 0;

			for (int j = 0; j < elemCountInRow; j++)
			{
				a_p[i] += A.val[A.rowPtr[i] + j] * p[A.colIndex[A.rowPtr[i] + j]];
			}

			alfa1 += r[i] * r_sop[i];
			alfa2 += a_p[i] * p_sop[i];
		}

		alfa = alfa1 / alfa2;

		// Пересчет x, r, r_sop, нахождение betta,
		// копирование в массивы predR и predR_sop

		double betta1 = 0;
		double betta2 = 0;

#pragma omp parallel for reduction(+:betta1, betta2)
		for (int i = 0; i < A.m; i++)
		{
			aT_p[i] = 0;
			for (int j = 0; j < A.colIndex.size(); j++)
			{
				if (A.colIndex[j] == i)
				{
					aT_p[i] += A.val[j] * p_sop[GetRowIndex(A, j)];
				}
			}

			if (i < n)
			{
				x[i] += alfa * p[i];
				r[i] -= alfa * a_p[i];
				r_sop[i] -= alfa * aT_p[i];

				betta1 += r[i] * r_sop[i];
				betta2 += predR[i] * predR_sop[i];

				predR[i] = r[i];
				predR_sop[i] = r_sop[i];
			}
		}

		if (A.m < n)
		{
#pragma omp parallel for reduction(+:betta1, betta2)
			for (int i = A.m; i < n; i++)
			{
				x[i] += alfa * p[i];
				r[i] -= alfa * a_p[i];
				r_sop[i] -= alfa * aT_p[i];

				betta1 += r[i] * r_sop[i];
				betta2 += predR[i] * predR_sop[i];

				predR[i] = r[i];
				predR_sop[i] = r_sop[i];
			}
		}

		betta = betta1 / betta2;

		// критерий останова

		if (count >= max_iter || abs(betta) < 0.000000001)
		{
			break;
		}

		double norma = 0;
#pragma omp parallel for reduction(+:norma)
		for (int i = 0; i < n; i++)
		{
			p[i] = r[i] + betta * p[i];
			p_sop[i] = r_sop[i] + betta * p_sop[i];

			norma += r[i] * r[i];
		}

		norma = sqrt(norma);
		if (norma < eps)
			break;

	} while (true);

	delete[] r, p, predR, r_sop, p_sop, predR_sop, aT_p, a_p;
}



