#include "stdafx.h"
#include "CppUnitTest.h"
#include "..\SlauSolver\SlauSolver.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace SlauSolverTests
{		
	TEST_CLASS(SluaSolverTests)
	{
	public:
		
		TEST_METHOD(IncrTest)
		{
			CSlauSolver solver;
			int res = solver.Incr(1);

			Assert::AreEqual(2, res);
		}

	};
}