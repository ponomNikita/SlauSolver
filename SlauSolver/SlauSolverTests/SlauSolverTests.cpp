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

		bool AreEqual(double *a, double *b, int n)
		{
			bool res = true;
			for (int i = 0; i < n; i++)
			{
				if (abs(a[i] - b[i]) > DOUBLE_EPS)
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

		TEST_METHOD(DiffVectorsTest)
		{
			double expected[] {-2, 0, 2};

			double a[] {1, 2, 3};
			double b[] {3, 2, 1};

			solver.Diff(a, b, 3);

			bool res = AreEqual(a, expected, 3);

			Assert::AreEqual(true, res);
		}

	};
}