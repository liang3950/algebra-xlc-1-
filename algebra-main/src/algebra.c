#include"algebra.h"
#include<stdio.h>
#include<math.h>

Matrix create_matrix(int row, int col){
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b){
    if(a.cols!=b.cols || a.rows!=b.rows){
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0,0);
    }else{
        Matrix c = create_matrix(a.rows, a.cols);
        for(int i=0;i<a.rows;i++){
            for(int j=0;j<a.cols;j++){
                c.data[i][j] = a.data[i][j] + b.data[i][j];
            }
        }
        return c;
    }
}

Matrix sub_matrix(Matrix a, Matrix b){
    if(a.cols!=b.rows || a.rows!=b.cols){
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0,0);
    }else{
        Matrix c = create_matrix(a.rows, a.cols);
        for(int i=0;i<a.rows;i++){
            for(int j=0;j<a.cols;j++){
                c.data[i][j] = a.data[i][j] - b.data[i][j];
            }
        }
        return c;
    }
}

Matrix mul_matrix(Matrix a, Matrix b){
    if(a.cols!=b.rows){
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0,0);
    }else{
        Matrix c = create_matrix(a.rows, b.cols);
        for(int i=0;i<a.rows;i++){
            for(int j=0;j<b.cols;j++){
                c.data[i][j] = 0;
                for(int k=0;k<a.cols;k++){
                    c.data[i][j] += a.data[i][k] * b.data[k][j];
                }
            }
        }
        return c;
    }
}

Matrix scale_matrix(Matrix a, double k){
    Matrix c = create_matrix(a.rows, a.cols);
    for(int i=0;i<a.rows;i++){
        for(int j=0;j<a.cols;j++){
            c.data[i][j] = a.data[i][j] * k;
        }
    }
    return c;
}

Matrix transpose_matrix(Matrix a){
    Matrix c = create_matrix(a.cols, a.rows);
    for(int i=0;i < a.rows; i++){
        for(int j=0;j < a.cols; j++){
            c.data[j][i] = a.data[i][j];
        }
    }
    return c;
}

double det_matrix(Matrix a){
    if(a.cols!=a.rows){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }else{
        if(a.rows==1){
            return a.data[0][0];
        }else{
            double det = 0;
            for(int i=0;i<a.cols;i++){
                Matrix c = create_matrix(a.rows-1, a.cols-1);
                for(int j=1;j<a.rows;j++){
                    for(int k=0;k<a.cols;k++){
                        if(k<i){
                            c.data[j-1][k] = a.data[j][k];
                        }else if(k>i){
                            c.data[j-1][k-1] = a.data[j][k];
                        }
                    }
                }
                det += a.data[0][i] * pow(-1, i) * det_matrix(c);
            }
            return det;
        }
    }
}

Matrix adjoint_matrix(Matrix a) {
    int n = a.rows;
    Matrix adj = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix sub = create_matrix(n - 1, n - 1);
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    if (k < i) {
                        if (l < j) {
                            sub.data[k][l] = a.data[k][l];
                        } else if (l > j) {
                            sub.data[k][l - 1] = a.data[k][l];
                        }
                    } else if (k > i) {
                        if (l < j) {
                            sub.data[k - 1][l] = a.data[k][l];
                        } else if (l > j) {
                            sub.data[k - 1][l - 1] = a.data[k][l];
                        }
                    }
                }
            }
            double det = det_matrix(sub);
            adj.data[j][i] = pow(-1, i + j) * det;
        }
    }
    return adj;
}

Matrix inv_matrix(Matrix a) {
    double det = det_matrix(a);
    if (det == 0) {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    } else {
        Matrix adj = adjoint_matrix(a);
        int n = a.rows;
        Matrix inv = create_matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inv.data[i][j] = adj.data[i][j] / det;
            }
        }
        return inv;
    }
}

int rank_matrix(Matrix a){
    int rank = 0;
    double eps = 1e-10;
    for (int i = 0; i < a.rows; ++i) {
        double maxEl = fabs(a.data[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < a.rows; ++k) {
            if (fabs(a.data[k][i]) > maxEl) {
                maxEl = fabs(a.data[k][i]);
                maxRow = k;
            }
        }
        if (maxEl < eps) {
            return rank;
        }
        if (maxRow != i) {
            double tmp;
            for(int k=0;k<a.cols;k++)
            {
                tmp = a.data[maxRow][k];
                a.data[maxRow][k] = a.data[i][k];
                a.data[i][k] = tmp;
            }
        }
        for (int k = i + 1; k < a.rows; ++k) {
            double c = -a.data[k][i] / a.data[i][i];
            for (int j = i; j < a.cols; ++j) {
                if (i == j) {
                    a.data[k][j] = 0;
                } else {
                    a.data[k][j] += c * a.data[i][j];
                }
            }
        }

        ++rank;
    }
    return rank;
}

double trace_matrix(Matrix a){
    if(a.rows!=a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }else{
        double trace = 0;
        for(int i=0;i<a.rows;i++){
            trace += a.data[i][i];
        }
        return trace;
    
    }
}

void print_matrix(Matrix a){
    for (int i = 0; i < a.rows; i++){
        for (int j = 0; j < a.cols; j++){
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}