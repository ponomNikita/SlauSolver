#include <vector>

using namespace std;

#ifdef SLAUSOLVER_EXPORTS
#define SLAUSOLVER_API __declspec(dllexport)
#else
#define SLAUSOLVER_API __declspec(dllimport)
#endif


struct CRSMatrix
{
	int n; // Число строк в матрице 
	int m; // Число столбцов в матрице 
	int nz; // Число ненулевых элементов в разреженной матрице 
	vector<double> val; // Массив значений матрицы по строкам 
	vector<int> colIndex; // Массив номеров столбцов 
	vector<int> rowPtr; // Массив индексов начала строк 
};

class SLAUSOLVER_API CSlauSolver {
public:
	CSlauSolver(void);

	void SLE_Solver_CRS_BICG(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count);

	double Dot(double *a, double *b, int n);

	// Сумма векторов с записью в первый вектор
	void Sum(double *a, double *b, int n);

	// Разность векторов с записью в первый вектор
	void Diff(double *a, double *b, int n);

	void Mult(CRSMatrix & A, double * b, double * res);
};

