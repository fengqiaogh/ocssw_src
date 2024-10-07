#include <iostream>

float** initializeMatrix(int n) {
    float** temp = new float*[n];
    for (int i = 0; i < n; i++)
        temp[i] = new float[n];
    return temp;
}

void input(float** M, int n) {
    std::cout << "Enter matrix: " << std::endl;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            std::cin >> M[i][j];
    std::cout << std::endl;
}

void printMatrix(float** M, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            std::cout << M[i][j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

float** add(float** M1, float** M2, int n) {
    float** temp = initializeMatrix(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            temp[i][j] = M1[i][j] + M2[i][j];
    return temp;
}

float** subtract(float** M1, float** M2, int n) {
    float** temp = initializeMatrix(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            temp[i][j] = M1[i][j] - M2[i][j];
    return temp;
}
