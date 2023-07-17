#include "s21_matrix.h"

// Создание матриц (create_matrix)
int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res = 0;
  if (result != NULL) {
    if (rows > 0 && columns > 0) {
      result->columns = columns;
      result->rows = rows;
      result->matrix = (double **)calloc(rows, sizeof(double *));
      if (result->matrix) {
        for (int i = 0; i < rows; i++) {
          result->matrix[i] = (double *)calloc(columns, sizeof(double));
          if (result->matrix[i] == NULL) {
            res = 1;
            break;
          }
        }
      }
    } else {
      matrix_to_NULL(result);
      res = 1;
    }
  } else {
    res = 1;
  }
  return res;
}

// Очистка матриц (remove_matrix)
void s21_remove_matrix(matrix_t *A) {
  if (A && A->matrix && A->rows > 0 && A->columns > 0) {
    for (int i = 0; i < A->rows && A->matrix[i] != NULL; i++) {
      free((double *)A->matrix[i]);
    }
    free((double **)A->matrix);
  }
  A->matrix = NULL;
  A->columns = 0;
  A->rows = 0;
}

// Сравнение матриц (eq_matrix)
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int result = 0;
  if (A->columns == B->columns && A->rows == B->rows) {
    for (int i = 0; i < A->rows; i++) {
      for (int y = 0; y < A->columns; y++) {
        if (floor(fabs((A->matrix[i][y] * pow(10, 7) -
                        B->matrix[i][y] * pow(10, 7)))) != 0.0) {
          result = 1;
          break;
        }
      }
    }
  } else {
    result = 1;
  }
  return !result;
}

// Сложение (sum_matrix)
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = 0;
  if (correct_matrix(A) && correct_matrix(B)) {
    if (A->columns == B->columns && A->rows == B->rows) {
      s21_create_matrix(A->rows, B->columns, result);
      for (int i = 0; i < result->rows; i++) {
        for (int y = 0; y < result->columns; y++) {
          result->matrix[i][y] = A->matrix[i][y] + B->matrix[i][y];
        }
      }
    } else
      error = 2;
  } else
    error = 1;
  return error;
}

// Вычитание матриц (sub_matrix)
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = 0;
  if (correct_matrix(A) && correct_matrix(B)) {
    if (A->columns == B->columns && A->rows == B->rows) {
      s21_create_matrix(A->rows, B->columns, result);
      for (int i = 0; i < result->rows; i++) {
        for (int y = 0; y < result->columns; y++) {
          result->matrix[i][y] = A->matrix[i][y] - B->matrix[i][y];
        }
      }
    } else
      error = 2;
  } else
    error = 1;
  return error;
}

// Умножение матрицы на число (mult_number).
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = 0;
  if (correct_matrix(A)) {
    matrix_t B;
    s21_create_matrix(A->rows, A->columns, &B);
    for (int i = 0; i < A->rows; i++) {
      for (int y = 0; y < A->columns; y++) {
        B.matrix[i][y] = number * A->matrix[i][y];
      }
    }
    s21_create_matrix(B.rows, B.columns, result);
    init_matrix(&B, result, B.rows + 1, B.columns + 1);
    s21_remove_matrix(&B);

  } else {
    error = 1;
  }
  return error;
}

// Умножение двух матриц (mult_matrix)
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = 0;
  if (!correct_matrix(A) || !correct_matrix(B)) {
    error = 1;
  } else if (A->columns != B->rows) {
    error = 2;
  } else {
    matrix_t C;
    int rows = (A->rows > B->rows) ? A->rows : B->rows;
    int columns = (A->columns > B->columns) ? A->columns : B->columns;
    s21_create_matrix(rows, columns, &C);
    for (int i = 0; i < C.rows; i++) {
      for (int y = 0; y < C.columns; y++) {
        for (int k = 0; k < B->rows; k++) {
          C.matrix[i][y] += A->matrix[i][k] * B->matrix[k][y];
        }
        C.matrix[i][y] = floor_matrix(C.matrix[i][y]);
      }
    }
    s21_create_matrix(A->rows, B->columns, result);
    init_matrix(&C, result, C.rows + 1, C.columns + 1);
    s21_remove_matrix(&C);
  }
  return error;
}

// Транспонирование матрицы (transpose)
int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = 0;
  if (correct_matrix(A)) {
    matrix_t B;
    int rows = A->rows, columns = A->columns;
    s21_create_matrix(columns, rows, &B);
    int sqr = (A->columns > A->rows) ? A->columns : A->rows;
    for (int i = 0; i < sqr; i++) {
      for (int y = 0; y < sqr; y++) {
        if (y == A->columns) break;
        B.matrix[y][i] = A->matrix[i][y];
      }
      if (i == A->rows - 1) break;
    }
    s21_create_matrix(columns, rows, result);
    init_matrix(&B, result, columns + 1, rows + 1);
    s21_remove_matrix(&B);
  } else {
    error = 1;
  }
  return error;
}

// Определитель матрицы (determinant)
int s21_determinant(matrix_t *A, double *result) {
  int error = 0;
  if (correct_matrix(A)) {
    if (A->rows == A->columns) {
      if (A->rows == 1) {
        *result = A->matrix[0][0];
      } else if (A->rows < 4) {
        *result = minor(A);
      } else {
        triangle(A, result);
      }
      *result = floor_matrix(*result);
      if (fabs(*result) < 1e-7) *result = 0;
    } else {
      error = 2;
      *result = 0;
    }
  } else {
    error = 1;
    *result = 0;
  }
  return error;
}

// Минор матрицы и матрица алгебраических дополнений (calc_complements)
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = 0;
  if (correct_matrix(A)) {
    if (A->rows == A->columns) {
      if (A->rows != 1) {
        matrix_t B, C;
        s21_create_matrix(A->rows - 1, A->columns - 1, &B);
        s21_create_matrix(A->rows, A->columns, &C);
        init_matrix(A, &C, A->rows + 1, A->columns + 1);
        for (int i = 0; i < A->rows; i++) {
          for (int y = 0; y < A->columns; y++) {
            double determinant = C.matrix[i][y];
            init_matrix(A, &B, i, y);
            s21_determinant(&B, &determinant);
            C.matrix[i][y] = pow(-1, (i + 1) + (y + 1)) * determinant;
          }
        }
        s21_create_matrix(A->rows, A->columns, result);
        init_matrix(&C, result, C.rows + 1, C.columns + 1);
        s21_remove_matrix(&C);
        s21_remove_matrix(&B);
      } else {
        error = 2;
      }
    } else {
      error = 2;
    }
  } else {
    error = 1;
  }
  return error;
}

// Обратная матрица (inverse_matrix)
int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = 0;
  if (correct_matrix(A)) {
    if (A->rows == A->columns && A->rows == 1) {
      s21_create_matrix(1, 1, result);
      result->matrix[0][0] = 1. / A->matrix[0][0];
    } else {
      double determinant = 0.0;
      s21_determinant(A, &determinant);
      if (determinant != 0) {
        matrix_t matrixminor, transposeMatrix;
        s21_calc_complements(A, &matrixminor);
        s21_transpose(&matrixminor, &transposeMatrix);
        s21_mult_number(&transposeMatrix, 1 / determinant, result);
        s21_remove_matrix(&matrixminor);
        s21_remove_matrix(&transposeMatrix);

      } else {
        error = 2;
      }
    }
  } else
    error = 1;
  return error;
}

//  доп функции
void print_matrix(matrix_t *matrix) {
  if (matrix->columns != 0 || matrix->rows != 0) {
    for (int i = 0; i < matrix->rows; i++) {
      for (int y = 0; y < matrix->columns; y++) {
        printf("\t%lf", matrix->matrix[i][y]);
      }
      printf("\n");
    }
  }
}

void matrix_to_NULL(matrix_t *A) {
  A->matrix = NULL;
  A->rows = 0;
  A->columns = 0;
}

int correct_matrix(matrix_t *A) {
  int res = 0;
  if (A == NULL || A->columns <= 0 || A->rows <= 0 || A->matrix == NULL) {
    res = 1;
  }
  return !res;
}

void init_matrix(matrix_t *A, matrix_t *result, int rows, int columns) {
  int k = 0, j = 0;
  for (int i = 0; i < A->rows; i++) {
    if (i == rows) i++;
    if (i == A->rows) break;
    for (int y = 0; y < A->columns; y++) {
      if (y == columns) y++;
      if (y == A->columns) break;
      result->matrix[k][j] = A->matrix[i][y];
      j++;
    }
    k++;
    j = 0;
  }
}

double minor(matrix_t *A) {
  double minor = 0;
  if (correct_matrix(A) && A->rows == A->columns) {
    if (A->rows == 2) {
      minor =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      double determinantSum = 1.0, determinantSub = 1.0;
      for (int k = 0; k < A->rows; k++) {
        for (int i = 0, y = 0; i < A->rows || y < A->columns; i++, y++) {
          if (y + k > A->rows - 1) {
            determinantSum *= A->matrix[i][(y + k) - A->rows];
          } else {
            determinantSum *= A->matrix[i][y + k];
          }
        }
        for (int i = 0, y = A->columns - 1; i < A->rows || y > 0; i++, y--) {
          if (y + k > A->rows - 1) {
            determinantSub *= A->matrix[i][(y + k) - A->rows];
          } else {
            determinantSub *= A->matrix[i][y + k];
          }
        }
        minor += determinantSum - determinantSub;
        determinantSum = 1.0;
        determinantSub = 1.0;
      }
    }
  }
  return minor;
}

void mult_string(matrix_t *A, double mul, int rows) {
  for (int i = 0; i < A->columns; i++) {
    A->matrix[rows][i] *= mul;
  }
}

void sub_string(matrix_t *A, int rows, int rows_sub) {
  for (int i = 0; i < A->columns; i++) {
    A->matrix[rows][i] = A->matrix[rows][i] - A->matrix[rows_sub][i];
  }
}

void div_string(matrix_t *A, double div, int rows) {
  for (int i = 0; i < A->columns; i++) {
    A->matrix[rows][i] /= div;
  }
}

void triangle(matrix_t *A, double *result) {
  matrix_t B;
  s21_create_matrix(A->rows, A->columns, &B);
  init_matrix(A, &B, A->rows + 1, A->columns + 1);
  for (int i = 0; i < A->rows; i++) {
    for (int y = i + 1; y < A->rows; y++) {
      if (B.matrix[y][i] == 0) continue;
      double mul = B.matrix[y][i] / B.matrix[i][i];
      mult_string(&B, mul, i);
      sub_string(&B, y, i);
      div_string(&B, mul, i);
    }
    *result = 1;
    for (int i = 0; i < B.rows; i++) {
      *result *= B.matrix[i][i];
    }
  }
  s21_remove_matrix(&B);
}

double floor_matrix(double A) {
  int sign = (A > 0) ? 1 : -1;
  if (A > 0) {
    A = sign * round(fabs(A) * pow(10, 7)) / pow(10, 7);
  } else {
    A = sign * round(fabs(A) * pow(10, 7)) / pow(10, 7);
  }
  return A;
}
