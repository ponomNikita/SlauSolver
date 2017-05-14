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
	tranposeA.m = A.n;
	tranposeA.n = A.m;
	tranposeA.val = vector<double>();
	tranposeA.colIndex = vector<int>();
	tranposeA.rowPtr = vector<int>();

	int k = 0;
	tranposeA.rowPtr.push_back(k);
	for (int i = 0; i < A.m; i++)
	{
		for (int j = 0; j < A.colIndex.size(); j++)
		{
			if (A.colIndex[j] == i)
			{
				tranposeA.val.push_back(A.val[j]);
				tranposeA.colIndex.push_back(GetRowIndex(A, j));
				k++;
			}
		}
		tranposeA.rowPtr.push_back(k);
	}
}

void CSlauSolver::Transpose1(CRSMatrix & A, CRSMatrix & At)
{
	At.n = A.m;
	At.m = A.n;
	At.nz = A.nz;
	At.val.resize(A.nz);
	At.colIndex.resize(A.nz);
	At.rowPtr.resize(A.m + 1);
	for (int i = 0; i < A.nz; i++)
		At.rowPtr[A.colIndex[i] + 1]++;

	int S = 0;
	int tmp_ = 0;

	for (int i = 1; i <= A.m; i++) {
		tmp_ = At.rowPtr[i];
		At.rowPtr[i] = S;
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
			int IIndex = At.rowPtr[RIndex + 1];
			At.val[IIndex] = v;
			At.colIndex[IIndex] = col;
			At.rowPtr[RIndex + 1]++;
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
	double* r = new double[n];
	double* p = new double[n];
	double* prevR = new double[n];
	double* prevP = new double[n];
	double* r_sop = new double[n];
	double* prevR_sop = new double[n];
	double* a_p = new double[n];
	double* aT_p = new double[n];
	double* p_sop = new double[n];
	double* prevP_sop = new double[n];

	CRSMatrix At;
	Transpose1(A, At);

#pragma omp parallel for
	for (int i = 0; i < n; i++) 
	{
		r[i] = b[i];
		r_sop[i] = r[i];
		x[i] = 0.0;
		p[i] = r[i];
		p_sop[i] = p[i];
	}

	for (; count < max_iter; count++) 
	{

#pragma omp parallel for
		for (int i = 0; i < n; i++) 
		{
			prevR[i] = r[i];
			prevP[i] = p[i];
			prevR_sop[i] = r_sop[i];
			prevP_sop[i] = p_sop[i];
		}

		Mult(A, prevP, a_p);
		Mult(At, prevP_sop, aT_p);

		double rr = Dot(prevR, prevR_sop, A.n);
		double App = Dot(a_p, prevP_sop, A.n);
		double alfa = rr / App;

#pragma omp parallel for
		for (int i = 0; i < n; i++) 
		{
			x[i] += alfa * p[i];
			r[i] = prevR[i] - alfa * a_p[i];
			r_sop[i] = prevR_sop[i] - alfa * aT_p[i]; 
		}

		rr = Dot(r, r_sop, A.n);
		double prevrr = Dot(prevR, prevR_sop, A.n);
		double beta = rr / prevrr;
		if (IsEnd(r, A.m, eps))
			break;

#pragma omp parallel for
		for (int i = 0; i < n; i++) 
		{
			p[i] = r[i] + beta * prevP[i];
			p_sop[i] = r_sop[i] + beta * prevP_sop[i];
		}
	}

	delete[] r, p, prevP, prevR, r_sop, prevR_sop, a_p, aT_p, p_sop, prevP_sop;
}



