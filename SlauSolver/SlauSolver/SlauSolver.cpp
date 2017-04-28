// SlauSolver.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "SlauSolver.h"


// This is an example of an exported variable
SLAUSOLVER_API int nSlauSolver=0;

// This is an example of an exported function.
SLAUSOLVER_API int fnSlauSolver(void)
{
	return 42;
}

// This is the constructor of a class that has been exported.
// see SlauSolver.h for the class definition
CSlauSolver::CSlauSolver()
{
	return;
}

int CSlauSolver::Incr(int a){
	return a + 1;
}
