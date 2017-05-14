#include "stdafx.h"
#include "SlauSolver.h"
#include <vector>

#include <omp.h>
using namespace std;

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

void CSlauSolver::Tranpose(CRSMatrix & A, CRSMatrix & tranposeA)
{
	tranposeA.n = A.m;
	tranposeA.m = A.n;
	tranposeA.nz = A.nz;
	tranposeA.val.resize(A.nz);
	tranposeA.colIndex.resize(A.nz);
	tranposeA.rowPtr.resize(A.m + 1);
	for (int i = 0; i < A.nz; i++)
		tranposeA.rowPtr[A.colIndex[i] + 1]++;

	int S = 0;
	int tmp_ = 0;

	for (int i = 1; i <= A.m; i++) {
		tmp_ = tranposeA.rowPtr[i];
		tranposeA.rowPtr[i] = S;
		S += tmp_;
	}

	for (int i = 0; i < A.n; i++) {
		int j1 = A.rowPtr[i];
		int j2 = A.rowPtr[i + 1];
		int col = i;
		for (int j = j1; j < j2; j++)
		{
			double v = A.val[j];
			int RIndex = A.colIndex[j];
			int IIndex = tranposeA.rowPtr[RIndex + 1];
			tranposeA.val[IIndex] = v;
			tranposeA.colIndex[IIndex] = col;
			tranposeA.rowPtr[RIndex + 1]++;
		}
	}
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
#pragma omp parallel for
	for (int i = 0; i < A.n; i++)
	{
		res[i] = 0;

		for (int j = A.rowPtr[i]; j < A.rowPtr[i + 1]; j++)
		{
			res[i] += A.val[j] * b[A.colIndex[j]];
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
	double * predX = new double[n]();

	double alfa, betta;

	CRSMatrix tranposeMatrix = {};

	Tranpose(A, tranposeMatrix);

	//GenerateSolution(x, n); // начальное приближение
	//SolveR(A, x, b, r, -1); // начальное r

#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
		r[i] = b[i];
		predR[i] = b[i];
		predR_sop[i] = b[i];
		r_sop[i] = b[i];
		p[i] = b[i];
		p_sop[i] = b[i];
		predX[i] = 0;
	}

	count = 0;

	for (; count < max_iter; count++)
	{
		double alfa1 = 0;
		double alfa2 = 0;

		// Вычисление alfa
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			a_p[i] = 0;

			for (int j = A.rowPtr[i]; j < A.rowPtr[i + 1]; j++)
			{
				a_p[i] += A.val[j] * p[A.colIndex[j]];
			}
		}

		for (int i = 0; i < n; i++)
		{
			alfa1 += r[i] * r_sop[i];
			alfa2 += a_p[i] * p_sop[i];
		}

		alfa = alfa1 / alfa2;

		// Пересчет x, r, r_sop, нахождение betta,
		// копирование в массивы predR и predR_sop

		double betta1 = 0;
		double betta2 = 0;


#pragma omp parallel for
			for (int i = 0; i < n; i++)
			{
				double tmp = 0;
				for (int j = tranposeMatrix.rowPtr[i]; j < tranposeMatrix.rowPtr[i + 1]; j++)
				{
					tmp += tranposeMatrix.val[j] * p_sop[tranposeMatrix.colIndex[j]];
				}

				x[i] += alfa * p[i];
				r[i] -= alfa * a_p[i];
				r_sop[i] -= alfa * tmp;
			}

			for (int i = 0; i < n; i++)
			{
				betta1 += r[i] * r_sop[i];
				betta2 += predR[i] * predR_sop[i];

				predR[i] = r[i];
				predR_sop[i] = r_sop[i];
			}

		betta = betta1 / betta2;

		// критерий останова

		if ( abs(betta) < 0.000000001)
		{
			break;
		}

		bool isEnd = true;
		bool endByXDiff = true;
		double norma = 0;
		for (int i = 0; i < n; i++)
		{
			p[i] = r[i] + betta * p[i];
			p_sop[i] = r_sop[i] + betta * p_sop[i];

			if (abs(x[i] - predX[i]) > eps)
			{
				endByXDiff = false;
			}
			norma += r[i] * r[i];
		}

		norma = sqrt(norma);

		isEnd = endByXDiff || norma < eps;

		if (isEnd)
			break;

	}

	delete[] r, p, predR, r_sop, p_sop, predR_sop, aT_p, a_p, predX;
}



