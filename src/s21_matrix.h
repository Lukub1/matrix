#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

// Основные
int s21_create_matrix(int rows, int columns, matrix_t *result);
// Очистка матриц (remove_matrix)
void s21_remove_matrix(matrix_t *A);
// Сравнение матриц (eq_matrix)
int s21_eq_matrix(matrix_t *A, matrix_t *B);
// Сложение (sum_matrix)
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// Вычитание матриц (sub_matrix)
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// Умножение матрицы на число (mult_number).
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
// Умножение двух матриц (mult_matrix)
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// Транспонирование матрицы (transpose)
int s21_transpose(matrix_t *A, matrix_t *result);
// Определитель матрицы (determinant)
int s21_determinant(matrix_t *A, double *result);
// Минор матрицы и матрица алгебраических дополнений (calc_complements)
int s21_calc_complements(matrix_t *A, matrix_t *result);
// Обратная матрица (inverse_matrix)
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

// доп функции
void print_matrix(matrix_t *matrix);
void matrix_to_NULL(matrix_t *A);
int correct_matrix(matrix_t *A);
void init_matrix(matrix_t *A, matrix_t *result, int rows, int columns);
double minor(matrix_t *A);
void mult_string(matrix_t *A, double mul, int rows);
void sub_string(matrix_t *A, int rows, int rows_sub);
void div_string(matrix_t *A, double div, int rows);
void triangle(matrix_t *A, double *result);
double floor_matrix(double A);

#endif  // SRC_S21_MATRIX_H_
