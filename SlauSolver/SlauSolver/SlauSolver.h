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
	void Sum(double *a, double *b, int n, double alfa = 1);
	
	// Сумма векторов с записью в другой вектор
	void Sum(double *a, double *b, double * res,  int n, double alfa = 1);

	// Разность векторов с записью в первый вектор
	void Diff(double *a, double *b, int n);

	void Mult(CRSMatrix & A, double * b, double * res);

	// Решение выражения вида Rk = Bk + alfa*A*Zk-1
	void SolveR(CRSMatrix & A, double * z, double * b, double * r, double alfa = 1);

	// Решение выражения вида Rk = Bk + alfa*A_transp*Zk-1
	void SolveRT(CRSMatrix & A, double * z, double * b, double * r, double alfa = 1);

	// Нахождение номера строки по индексу элемента в массиве val
	int GetRowIndex(CRSMatrix & A, int index);

private: 
	void GenerateSolution(double * x, int n);
	void Copy(double * a, double * copyA, int n);
	bool IsEnd(double * x0, double * x, int n, double eps);
};

