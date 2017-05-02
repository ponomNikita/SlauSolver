#include "stdafx.h"
#include "CppUnitTest.h"
#include "..\SlauSolver\SlauSolver.h"

#define DOUBLE_EPS 0.000000001

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace SlauSolverTests
{		
	TEST_CLASS(SluaSolverTests)
	{
	private:
		CSlauSolver solver;

		bool AreEqual(double *a, double *b, int n, double eps = DOUBLE_EPS)
		{
			bool res = true;
			for (int i = 0; i < n; i++)
			{
				if (abs(a[i] - b[i]) > eps)
				{
					res = false;
					break;
				}
			}

			return res;
		}

	public:
		
		TEST_METHOD(DotTest)
		{
			double expected = 10;

			double a[] {1, 2, 3};
			double b[] {3, 2, 1};

			double res = solver.Dot(a, b, 3);

			Assert::AreEqual(expected, res);
		}

		TEST_METHOD(SumVectorsTest)
		{
			double expected[] {4, 4, 4};

			double a[] {1, 2, 3};
			double b[] {3, 2, 1};

			solver.Sum(a, b, 3);

			bool res = AreEqual(a, expected, 3);

			Assert::AreEqual(true, res);
		}

		TEST_METHOD(SumVectorsTest2)
		{
			double expected[] {7, 6, 5};

			double a[] {1, 2, 3};
			double b[] {3, 2, 1};

			double * resVector = new double[3];

			solver.Sum(a, b, resVector, 3, 2);

			bool res = AreEqual(resVector, expected, 3);

			Assert::AreEqual(true, res);
		}

		TEST_METHOD(DiffVectorsTest)
		{
			double expected[] {-2, 0, 2};

			double a[] {1, 2, 3};
			double b[] {3, 2, 1};

			solver.Diff(a, b, 3);

			bool res = AreEqual(a, expected, 3);

			Assert::AreEqual(true, res);
		}

		TEST_METHOD(MultCRSMatrixOnVectorTest)
		{
			/*
				10, 0, 0, -2, 0,
				3, 9, 0, 0, 0, 3,
				0, 7, 8, 7, 0, 0,
				3, 0, 8, 7, 5, 0,
				0, 8, 0, 9, 9, 13,
				0, 4, 0, 0, 2, -1
			*/
			vector<double> val = 
			{
				10, -2, 3, 9, 3, 7, 8, 7, 3, 8, 7, 5, 8, 9, 9, 13, 4, 2, -1
			}; 
			vector<int> colIndex = 
			{
				0, 4, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
			};
			vector<int> rowPtr = 
			{
				0, 2, 5, 8, 12, 16, 19
			};

			CRSMatrix matrix = {};
			matrix.n = 6;
			matrix.m = 6;
			matrix.val = val;
			matrix.rowPtr = rowPtr;
			matrix.colIndex = colIndex;

			double b[] = 
			{
				1, 1, 1, 1, 1, 1
			};

			double expectedResult[] = 
			{
				8, 15, 22, 23, 39, 5
			};

			double * res = new double[6];

			solver.Mult(matrix, b, res);

			Assert::IsTrue(AreEqual(res, expectedResult, 6));
		
		}

		TEST_METHOD(MultCRSMatrixWithZeroRowOnVectorTest)
		{
			/*
			10, 0, 0, -2, 0,
			0, 0, 0, 0, 0, 0,
			0, 7, 8, 7, 0, 0,
			3, 0, 8, 7, 5, 0,
			0, 8, 0, 9, 9, 13,
			0, 4, 0, 0, 2, -1
			*/
			vector<double> val =
			{
				10, -2, 7, 8, 7, 3, 8, 7, 5, 8, 9, 9, 13, 4, 2, -1
			};
			vector<int> colIndex =
			{
				0, 4, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
			};
			vector<int> rowPtr =
			{
				0, 2, 2, 5, 9, 13, 16
			};

			CRSMatrix matrix = {};
			matrix.n = 6;
			matrix.m = 6;
			matrix.val = val;
			matrix.rowPtr = rowPtr;
			matrix.colIndex = colIndex;

			double b[] =
			{
				1, 1, 1, 1, 1, 1
			};

			double expectedResult[] =
			{
				8, 0, 22, 23, 39, 5
			};

			double * res = new double[6];

			solver.Mult(matrix, b, res);

			Assert::IsTrue(AreEqual(res, expectedResult, 6));

		}

		TEST_METHOD(SolveRTest)
		{
			/*
			10, 0, 0, -2, 0,
			3, 9, 0, 0, 0, 3,
			0, 7, 8, 7, 0, 0,
			3, 0, 8, 7, 5, 0,
			0, 8, 0, 9, 9, 13,
			0, 4, 0, 0, 2, -1
			*/
			vector<double> val =
			{
				10, -2, 3, 9, 3, 7, 8, 7, 3, 8, 7, 5, 8, 9, 9, 13, 4, 2, -1
			};
			vector<int> colIndex =
			{
				0, 4, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
			};
			vector<int> rowPtr =
			{
				0, 2, 5, 8, 12, 16, 19
			};

			CRSMatrix matrix = {};
			matrix.n = 6;
			matrix.m = 6;
			matrix.val = val;
			matrix.rowPtr = rowPtr;
			matrix.colIndex = colIndex;

			double b[] =
			{
				1, 1, 1, 1, 1, 1
			};

			double expectedResult[] =
			{
				-7, -14, -21, -22, -38, -4
			};

			double alfa = -1;

			double r[] =
			{
				1, 1, 1, 1, 1, 1
			};

			solver.SolveR(matrix, b, r, alfa);

			Assert::IsTrue(AreEqual(r, expectedResult, 6));

		}

	};
}