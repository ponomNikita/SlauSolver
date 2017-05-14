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

		bool isInfinite(const double &value)
		{
			double max_value = numeric_limits<double>::max();
			double min_value = -max_value;

			return !(min_value <= value && value <= max_value);
		}

		bool AreEqual(double *a, double *b, int n, double eps = DOUBLE_EPS)
		{
			bool res = true;
			for (int i = 0; i < n; i++)
			{
				if (abs(a[i] - b[i]) > eps || isInfinite(a[i]) || isInfinite(b[i]))
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
				10, 0, 0, -2, 0, 0
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
				0, 3, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
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

		TEST_METHOD(MultTransposeCRSMatrixOnVectorTest)
		{
			/*
			10, 0, 0, -2, 0, 0
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
				0, 3, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
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
				16, 28, 16, 21, 16, 15
			};

			double * res = new double[6];

			solver.MultTransp(matrix, b, res);

			Assert::IsTrue(AreEqual(res, expectedResult, 6));

		}

		TEST_METHOD(MultCRSMatrixWithZeroRowOnVectorTest)
		{
			/*
			10, 0, 0, -2, 0, 0
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
				0, 3, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
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
			10, 0, 0, -2, 0, 0
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
				0, 3, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
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

			solver.SolveR(matrix, b, r, r, alfa);

			Assert::IsTrue(AreEqual(r, expectedResult, 6));

		}

		TEST_METHOD(SolveRTranspTest)
		{
			/*
			10, 0, 0, -2, 0, 0
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
				0, 3, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
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
				//-7, -14, -21, -22, -38, -4
				-15, -27, -15, -20, -15, -14
			};

			double alfa = -1;

			double r[] =
			{
				1, 1, 1, 1, 1, 1
			};

			solver.SolveRT(matrix, b, r, r, alfa);

			Assert::IsTrue(AreEqual(r, expectedResult, 6));

		}

		TEST_METHOD(SLE_Solver_CRS_BICG_Test)
		{
			/*
			2, 1, 1
			1, -1, 0,
			3, -1, 2			
			*/

			vector<double> val =
			{
				2, 1,  1, 1, -1, 3, -1, 2
			};
			vector<int> colIndex =
			{
				0, 1, 2, 0, 1, 0, 1, 2
			};
			vector<int> rowPtr =
			{
				0, 3, 5, 8 
			};

			CRSMatrix matrix = {};
			matrix.n = 3;
			matrix.m = 3;
			matrix.nz = val.size();
			matrix.val = val;
			matrix.rowPtr = rowPtr;
			matrix.colIndex = colIndex;

			double b[] = {
				2, -2, 2
			};

			double expectedX[] = {
				-1, 1, 3
			};

			double * x = new double[3];

			int count = 0;

			solver.SLE_Solver_CRS_BICG(matrix, b, 0.001, 100, x, count);

			Assert::IsTrue(AreEqual(x, expectedX, 3, 0.001));
		}

		TEST_METHOD(GetRowIndexTest)
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
			matrix.nz = val.size();
			matrix.val = val;
			matrix.rowPtr = rowPtr;
			matrix.colIndex = colIndex;

			int res1 = solver.GetRowIndex(matrix, 0);
			int res2 = solver.GetRowIndex(matrix, 2);
			int res3 = solver.GetRowIndex(matrix, 9);
			int res4 = solver.GetRowIndex(matrix, 14);

			int expectedRes1 = 0;
			int expectedRes2 = 1;
			int expectedRes3 = 3;
			int expectedRes4 = 4;

			Assert::AreEqual(expectedRes1, res1);
			Assert::AreEqual(expectedRes2, res2);
			Assert::AreEqual(expectedRes3, res3);
			Assert::AreEqual(expectedRes4, res4);

		}

		TEST_METHOD(TransposeMatrixTest)
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
				0, 3, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5
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

			double valT[] =
			{
				10, 3, 3, 9, 7, 8, 4, 8, 8, -2, 7, 7, 9, 5, 9, 2, 3, 13, -1
			};
			int colIndexT[] =
			{
				0, 1, 3, 1, 2, 4, 5, 2, 3, 0, 2, 3, 4, 3, 4, 5, 1, 4, 5
			};
			int rowPtrT[] =
			{
				0, 3, 7, 9, 13, 16, 19
			};

			double * expectedVal = new double[19];
			int * expectedColIndex = new int[19];
			int * expectedrowPtr = new int[7];

			CRSMatrix transposeMatrix = {};

			solver.Tranpose(matrix, transposeMatrix);

			for (int i = 0; i < 19; i++)
			{
				expectedVal[i] = transposeMatrix.val[i];
			}

			Assert::IsTrue(AreEqual(valT, expectedVal, 19, 0.0001));
		}

		TEST_METHOD(SLE_Solver_CRS_BICG_Test2)
		{
			/*
			2, 3, -4, 1
			1, 0, 0, 1
			-1, -1, 2, -1
			2, -1, -1, 2
			*/

			vector<double> val =
			{
				2, 3, -4, 1, 1, 1, -1, -1, 2, -1, 2, -1, -1, 2
			};
			vector<int> colIndex =
			{
				0, 1, 2, 3, 0, 3, 0, 1, 2, 3, 0, 1, 2, 3
			};
			vector<int> rowPtr =
			{
				0, 4, 6, 10, 14
			};

			CRSMatrix matrix = {};
			matrix.n = 4;
			matrix.m = 4;
			matrix.nz = val.size();
			matrix.val = val;
			matrix.rowPtr = rowPtr;
			matrix.colIndex = colIndex;

			double b[] = {
				0, -1, 1, -2
			};

			double expectedX[] = {
				1, 0, 0, -2
			};

			double * x = new double[4];

			int count = 0;

			solver.SLE_Solver_CRS_BICG(matrix, b, 0.001, 100, x, count);

			Assert::IsTrue(AreEqual(x, expectedX, 4, 0.001));
		}

		TEST_METHOD(SolveRTest2)
		{
			/*
			2, 3, -4, 1
			1, 0, 0, 4
			-1, -1, 2, -1
			2, -1, -1, 2
			*/

			vector<double> val =
			{
				2, 3, -4, 1, 1, 4, -1, -1, 2, -1, 2, -1, -1, 2
			};
			vector<int> colIndex =
			{
				0, 1, 2, 3, 0, 3, 0, 1, 2, 3, 0, 1, 2, 3
			};
			vector<int> rowPtr =
			{
				0, 4, 6, 10, 14
			};

			CRSMatrix matrix = {};
			matrix.n = 4;
			matrix.m = 4;
			matrix.val = val;
			matrix.rowPtr = rowPtr;
			matrix.colIndex = colIndex;

			double b[] =
			{
				1, 1, 1, 1
			};

			double expectedResult[] =
			{
				-1, -4, 2, -1
			};

			double alfa = -1;

			double r[] =
			{
				1, 1, 1, 1
			};

			solver.SolveR(matrix, b, r, r, alfa);

			Assert::IsTrue(AreEqual(r, expectedResult, 4));

		}

		TEST_METHOD(SolveRTranspTest2)
		{
			/*
			2, 3, -4, 1
			1, 0, 0, 4
			-1, -1, 2, -1
			2, -1, -1, 2
			*/

			vector<double> val =
			{
				2, 3, -4, 1, 1, 4, -1, -1, 2, -1, 2, -1, -1, 2
			};
			vector<int> colIndex =
			{
				0, 1, 2, 3, 0, 3, 0, 1, 2, 3, 0, 1, 2, 3
			};
			vector<int> rowPtr =
			{
				0, 4, 6, 10, 14
			};

			CRSMatrix matrix = {};
			matrix.n = 4;
			matrix.m = 4;
			matrix.val = val;
			matrix.rowPtr = rowPtr;
			matrix.colIndex = colIndex;

			double b[] =
			{
				1, 1, 1, 1
			};

			double expectedResult[] =
			{
				-3, 0, 4, -5
			};

			double alfa = -1;

			double r[] =
			{
				1, 1, 1, 1
			};

			solver.SolveRT(matrix, b, r, r, alfa);

			Assert::IsTrue(AreEqual(r, expectedResult, 4));

		}

	};
}