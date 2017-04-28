#include "stdafx.h"
#include "SlauSolver.h"

CSlauSolver::CSlauSolver()
{
	return;
}

void CSlauSolver::SLE_Solver_CRS_BICG(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{

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

void CSlauSolver::Sum(double *a, double *b, int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] += b[i];
	}
}

void CSlauSolver::Diff(double *a, double *b, int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] -= b[i];
	}
}

