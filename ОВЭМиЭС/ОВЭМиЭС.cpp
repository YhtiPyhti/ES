#include <iostream>
#include <ctime>
#include <immintrin.h>


using namespace std;

int checkInt(string str);

void initMatrix(float* matrix[], int M, int N);

int DGEMM(float* matrix2[], int M2, int N2, float* matrix1[], int M1, int N1);
int DGEMM_opt_1(float* matrix2[], int M2, int N2, float* matrix1[], int M1, int N1);
int DGEMM_opt_2(float* matrix2[], int M2, int N2, float* matrix1[], int M1, int N1, int block_size);
int DGEMM_opt_3(float* matrix2[], int M2, int N2, float* matrix1[], int M1, int N1);

int main(int argc, char* argv[])
{


    int LenM1_Str = checkInt(argv[1]), LenM1_Stl = checkInt(argv[2]);
    int LenM2_Str = checkInt(argv[3]), LenM2_Stl = checkInt(argv[4]), block = checkInt(argv[5]);

    if (LenM1_Str < 0 || LenM2_Str < 0 || LenM1_Stl < 0 || LenM2_Stl < 0 || block < 0) return -1;
    srand(time(0));

    cout << "Raw Matrix 1: " << LenM1_Str << " Column Matrix 1: " << LenM1_Stl << endl;
    cout << "Raw Matrix 2: " << LenM2_Str << " Column Matrix 2: " << LenM2_Stl << endl;

    float** matrix1 = new float* [LenM1_Str];
    for (size_t count = 0; count < LenM1_Str; count++)
        matrix1[count] = new float[LenM1_Stl];

    float** matrix2 = new float* [LenM2_Str];
    for (size_t count = 0; count < LenM2_Str; count++)
        matrix2[count] = new float[LenM2_Stl];


    initMatrix(matrix1, LenM1_Str, LenM1_Stl);
    initMatrix(matrix2, LenM2_Str, LenM2_Stl);
    cout << "Matrix init" << endl;

    DGEMM(matrix2, LenM2_Str, LenM2_Stl, matrix1, LenM1_Str, LenM1_Stl);
    DGEMM_opt_1(matrix2, LenM2_Str, LenM2_Stl, matrix1, LenM1_Str, LenM1_Stl);
    DGEMM_opt_2(matrix2, LenM2_Str, LenM2_Stl, matrix1, LenM1_Str, LenM1_Stl, block);
    DGEMM_opt_3(matrix2, LenM2_Str, LenM2_Stl, matrix1, LenM1_Str, LenM1_Stl);

    for (size_t count = 0; count < LenM1_Str; count++)
        delete[] matrix1[count];

    for (size_t count = 0; count < LenM1_Str; count++)
        delete[] matrix2[count];


    return 0;
}

int checkInt(string str) {
    int number = 0;
    for (size_t i = 0, lenght = str.length(); i < lenght; i++)
    {
        if (str[i] < 48 || str[i] > 58) return -1;
        number *= 10;
        number += str[i] - 48;
    }
    return number;
}

void initMatrix(float* matrix[], int M, int N) {
    for (size_t count_row = 0; count_row < M; count_row++)
        for (size_t count_column = 0; count_column < N; count_column++)
            matrix[count_row][count_column] = (rand() % 10 + 1) / float((rand() % 10 + 1));
}

int DGEMM(float* matrix2[], int M2, int N2, float* matrix1[], int M1, int N1) {

    if (N2 != M1) return -1;

    float** matrix_result = new float* [M1];
    for (size_t count = 0; count < M1; count++)
        matrix_result[count] = new float[N2];
    cout << "Matrix create" << endl;
    for (size_t count_row = 0; count_row < M1; count_row++)
        for (int count_column = 0; count_column < N2; count_column++)
            matrix_result[count_row][count_column] = 0;

    clock_t start = clock();

    int i = 0;

    for (size_t count_row = 0; count_row < M1; count_row++) {
        for (size_t count_column = 0; count_column < N2; count_column++) {
            i = 0;
            while (i < N1)
            {
                matrix_result[count_row][count_column] += matrix1[count_row][i] * matrix2[i][count_column];
                i++;
            }
        }
    }

    clock_t end = clock();

    double sec = (double)(end - start) / CLOCKS_PER_SEC;

    cout << "Time no opt: " << sec << endl;

    for (size_t count = 0; count < M1; count++)
        delete[] matrix_result[count];

    return -1;
}

int DGEMM_opt_1(float* matrix2[], int M2, int N2, float* matrix1[], int M1, int N1) {

    if (N2 != M1) return -1;

    float** matrix_result = new float* [M1];
    for (size_t count = 0; count < M1; count++)
        matrix_result[count] = new float[N2];

    for (size_t count_row = 0; count_row < M1; count_row++)
        for (size_t count_column = 0; count_column < N2; count_column++)
            matrix_result[count_row][count_column] = 0;

    clock_t start = clock();

    int i = 0;

    for (size_t count_row = 0; count_row < M1; count_row++) {
        for (size_t count_column = 0; count_column < N2; count_column++) {
            for (size_t i = 0; i < N1; i++) {
                matrix_result[count_row][count_column] += matrix1[count_row][count_column] * matrix2[count_column][i];
            }
        }
    }

    clock_t end = clock();

    double sec = (double)(end - start) / CLOCKS_PER_SEC;

    cout << "Time column: " << sec << endl;

    for (int count = 0; count < M1; count++)
        delete[] matrix_result[count];

    return -1;
}

int DGEMM_opt_2(float* matrix2[], int M2, int N2, float* matrix1[], int M1, int N1, int block_size) {

    if (N2 != M1) return -1;
    float sum = 0;

    float** matrix_result = new float* [M1];
    for (size_t count = 0; count < M1; count++)
        matrix_result[count] = new float[N2];

    for (size_t count_row = 0; count_row < M1; count_row++)
        for (size_t count_column = 0; count_column < N2; count_column++)
            matrix_result[count_row][count_column] = 0;

    clock_t start = clock();

    int i = 0, lenght = 0;

    for (size_t count_row = 0; count_row < M1; count_row += block_size) {
        for (size_t count_column = 0; count_column < N2; count_column += block_size) {
            for (size_t i = 0; i < N1; i++) {
                lenght = (count_row + block_size) > N2 ? N2 : count_row + block_size;
                for (size_t j = count_row; j < lenght; j++) {
                    sum = 0;
                    lenght = (count_column + block_size) > N2 ? N2 : count_column + block_size;
                    for (size_t k = count_column; k < lenght; k++)
                        sum += matrix1[i][k] * matrix1[k][j];

                    matrix_result[i][j] += sum;
                }
            }
        }
    }

    clock_t end = clock();

    double sec = (double)(end - start) / CLOCKS_PER_SEC;

    cout << "Time block: " << sec << endl;

    for (int count = 0; count < M1; count++)
        delete[] matrix_result[count];

    return -1;
}

int DGEMM_opt_3(float* matrix2[], int M2, int N2, float* matrix1[], int M1, int N1) {

    if (N2 != M1) return -1;

    float** matrix_result = new float* [M1];
    for (size_t count = 0; count < M1; count++)
        matrix_result[count] = new float[N2];

    for (size_t count_row = 0; count_row < M1; count_row++)
        for (int count_column = 0; count_column < N2; count_column++)
            matrix_result[count_row][count_column] = 0;

    clock_t start = clock();

    int i = 0;

    for (int count_row = 0; count_row < M1; count_row++)
    {
        float* c = matrix_result[count_row];
        for (int j = 0; j < N2; j += 8)
            _mm256_storeu_ps(c + j, _mm256_setzero_ps());
        for (int k = 0; k < N1; ++k)
        {
            const float* b = matrix2[k];
            __m256 a = _mm256_set1_ps(matrix1[i][k]);
            for (int j = 0; j < N1; j += 16)
            {
                _mm256_storeu_ps(c + j + 0, _mm256_fmadd_ps(a,
                    _mm256_loadu_ps(b + j + 0), _mm256_loadu_ps(c + j + 0)));
                _mm256_storeu_ps(c + j + 8, _mm256_fmadd_ps(a,
                    _mm256_loadu_ps(b + j + 8), _mm256_loadu_ps(c + j + 8)));
            }
        }
    }


    clock_t end = clock();

    double sec = (double)(end - start) / CLOCKS_PER_SEC;

    cout << "Time vector: " << sec << endl;
    for (size_t count = 0; count < M1; count++)
        delete[] matrix_result[count];

    return -1;
}