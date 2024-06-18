#include <iostream>
#include <vector>
#include <cmath>

#define GRID_WIDTH 11

typedef std::vector<double> RealVector;
typedef std::vector<RealVector> RealMatrix;

// 计算向量p范数
double norm(const RealVector &x, int p) {
	int len = x.size();
	double ans = 0;
	for (int i = 0; i < len; i++)
		ans += pow(abs(x[i]), p);
	return pow(ans, 1 / p);
}

// 计算向量x-y的p范数
double norm_sub(const RealVector &x, const RealVector &y, int p) {
	int len = x.size();
	double ans = 0;
	for (int i = 0; i < len; i++)
		ans += pow(abs(x[i] - y[i]), p);
	return pow(ans, 1 / p);
}

// Jacobi迭代法
void solve_jacobi(RealMatrix &a, RealVector &x, RealVector &b, int maxn, double eps) {
    int len = x.size();
	for (int k = 0; k < maxn; k++) {
		RealVector backup(x.begin(), x.end());
		for (int i = 0; i < len; i++) {
			double sum = 0;
			for (int j = 0; j < len; j++)
				sum += a[i][j] * backup[j];
			x[i] = backup[i] + 1 / a[i][i] * (b[i] - sum);
		}
		
		printf("iterator = %d, error = %e\n", k + 1, norm_sub(x, backup, 1));
		// printf("x = ");
		
		for (int i = 0; i < len; i++) {
			printf("+");
			for (int j = 0; j < GRID_WIDTH; j++)
				printf("-");
		}
		printf("+\n");
		for (int i = 0; i < len; i++)
			printf("|%11.4e", x[i]); // GRID_WIDTH
		printf("|\n");
		for (int i = 0; i < len; i++) {
			printf("+");
			for (int j = 0; j < GRID_WIDTH; j++)
				printf("-");
		}
		printf("+\n\n");

		if (norm_sub(x, backup, 1) < eps)
			break;
	}
	return ;
}

// sor迭代法，omega = 1时为GS迭代法
void solve_sor(RealMatrix &a, RealVector &x, RealVector &b, int maxn, double omega, double eps) {
    int len = x.size();
	for (int k = 0; k < maxn; k++) {
		RealVector backup(x.begin(), x.end());
		for (int i = 0; i < len; i++) {
			double lhs = 0, rhs = 0;
			for (int j = 0; j < i; j++)
				lhs += a[i][j] * x[j];
			for (int j = i; j < len; j++)
				rhs += a[i][j] * x[j];
			x[i] = x[i] + omega / a[i][i] * (b[i] - lhs - rhs);
		}

		printf("iterator = %d, error = %e\n", k + 1, norm_sub(x, backup, 1));
		// printf("x = ");
		
		for (int i = 0; i < len; i++) {
			printf("+");
			for (int j = 0; j < GRID_WIDTH; j++)
				printf("-");
		}
		printf("+\n");
		for (int i = 0; i < len; i++)
			printf("|%11.4e", x[i]); // GRID_WIDTH
		printf("|\n");
		for (int i = 0; i < len; i++) {
			printf("+");
			for (int j = 0; j < GRID_WIDTH; j++)
				printf("-");
		}
		printf("+\n\n");
		
		if (norm_sub(x, backup, 1) < eps)
			break;
	}
	return ;
}

int main()
{
    int n;
	printf("Please enter a number for matrix size: ");
    scanf("%d", &n);

    RealMatrix arr(n, RealVector(n, 0));
	printf("Please enter values for matrix A: \n");
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < n; j++) 
            scanf("%lf", &arr[i][j]);
    
    RealVector bv(n, 0);
	printf("Please enter values for vector b: \n");
    for (int i = 0; i < n; i++)  
        scanf("%lf", &bv[i]);
    
    RealVector xv(n, 0);
	printf("Please enter values for initial vector x: \n");
	for (int i = 0; i < n; i++) 
    	scanf("%lf", &xv[i]); 

	int maxn;
	printf("Please enter a number for maximum time of iteration: ");
	scanf("%d", &maxn);

	double epsilon;
	printf("Please enter a number for maximum error: ");
	scanf("%lf", &epsilon);

	int mode = 1;
	printf("If you enter 0, the program will use Jacobi iteration method, ");
	printf("otherwise SOR iteration method will be used. \n");
	printf("Please choose a method: ");
	scanf("%d", &mode);

	if (!mode) {
		solve_jacobi(arr, xv, bv, maxn, epsilon);
	}
	else {
		double w;
		printf("Please enter SOR factor(0 < w < 2): ");
		scanf("%lf", &w);
		solve_sor(arr, xv, bv, maxn, w, epsilon);
	}
	return 0;
}