#include <vector>

using namespace std;

#ifdef SLAUSOLVER_EXPORTS
#define SLAUSOLVER_API __declspec(dllexport)
#else
#define SLAUSOLVER_API __declspec(dllimport)
#endif


struct CRSMatrix
{
	int n; // ����� ����� � ������� 
	int m; // ����� �������� � ������� 
	int nz; // ����� ��������� ��������� � ����������� ������� 
	vector<double> val; // ������ �������� ������� �� ������� 
	vector<int> colIndex; // ������ ������� �������� 
	vector<int> rowPtr; // ������ �������� ������ ����� 
};

class SLAUSOLVER_API CSlauSolver {
public:
	CSlauSolver(void);

	void SLE_Solver_CRS_BICG(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count);

	double Dot(double *a, double *b, int n);

	// ����� �������� � ������� � ������ ������
	void Sum(double *a, double *b, int n);

	// �������� �������� � ������� � ������ ������
	void Diff(double *a, double *b, int n);

	void Mult(CRSMatrix & A, double * b, double * res);
};

