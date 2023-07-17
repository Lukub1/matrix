#include "s21_matrix.h"

#include <check.h>

void filling_matrix(matrix_t *A, double number);
void filling_diag(matrix_t *A);
void filling_arr(matrix_t *A, double *number);

START_TEST(test_s21_create_matrix) {
  // не корректные матрицы
  matrix_t A;
  int res = s21_create_matrix(0, 5, &A);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);
  res = s21_create_matrix(5, 0, &A);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);
  res = s21_create_matrix(2, -2, &A);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);
  res = s21_create_matrix(-2, 5, &A);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);

  // корректные матрицы
  res = s21_create_matrix(5, 5, &A);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&A);
  res = s21_create_matrix(1, 5, &A);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&A);
  res = s21_create_matrix(1, 1, &A);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_s21_eq_matrix) {
  matrix_t A, B, C, D, E;

  s21_create_matrix(4, 3, &A);
  s21_create_matrix(4, 3, &B);
  s21_create_matrix(3, 3, &C);
  s21_create_matrix(3, 3, &D);
  s21_create_matrix(-3, 3, &E);

  filling_matrix(&A, 0.1234567);
  filling_matrix(&B, 0.1234568);
  filling_matrix(&C, 0.12345678);
  filling_matrix(&D, 0.12345679);

  //  матрицы с разными значениями и точностью 7
  int res = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res, 0);

  //  матрицы разных размеров
  res = s21_eq_matrix(&A, &C);
  ck_assert_int_eq(res, 0);

  //  матрицы с одинаовыми значенями
  res = s21_eq_matrix(&A, &A);
  ck_assert_int_eq(res, 1);

  //  матрицы с разными значениями и точностью больше 7
  res = s21_eq_matrix(&B, &D);
  ck_assert_int_eq(res, 0);

  //  корректная и некорректная матрица
  res = s21_eq_matrix(&A, &E);
  ck_assert_int_eq(res, 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&D);
  s21_remove_matrix(&E);
}
END_TEST

START_TEST(test_s21_sum_matrix) {
  matrix_t A, B, C, result, resultoriginal;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 3, &B);
  double numberA[] = {1,   2, 3, 4.8, 6,   8, 5.5, 7, -1, 2, 3,   4, -6,
                      8.5, 5, 7, -6,  8.4, 5, 7,   1, 0,  5, 1.2, 3};
  filling_arr(&A, numberA);
  double numberB[] = {1, 2, -3, 8, -5.5, 7, 3, -4, 6, 7, 6, -8.4};
  filling_arr(&B, numberB);
  // корректные матрицы одного размера
  int res = s21_sum_matrix(&A, &A, &result);
  ck_assert_int_eq(res, 0);
  s21_create_matrix(5, 5, &resultoriginal);
  double numberoriginalA[] = {2,  4,  6, 9.6, 12, 16,  11, 14,  -2,
                              4,  6,  8, -12, 17, 10,  14, -12, 16.8,
                              10, 14, 2, 0,   10, 2.4, 6};
  filling_arr(&resultoriginal, numberoriginalA);
  res = s21_eq_matrix(&result, &resultoriginal);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&resultoriginal);

  res = s21_sum_matrix(&B, &B, &result);
  ck_assert_int_eq(res, 0);
  s21_create_matrix(4, 3, &resultoriginal);
  double numberoriginalB[] = {2, 4, -6, 16, -11, 14, 6, -8, 12, 14, 12, -16.8};
  filling_arr(&resultoriginal, numberoriginalB);
  res = s21_eq_matrix(&result, &resultoriginal);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&resultoriginal);

  // корректные матрицы разных размеров
  res = s21_sum_matrix(&A, &B, &result);
  ck_assert_uint_eq(res, 2);
  s21_remove_matrix(&result);

  // некорректная матрица
  s21_create_matrix(0, -3, &C);
  res = s21_sum_matrix(&C, &C, &result);
  ck_assert_uint_eq(res, 1);
  res = s21_sum_matrix(&C, &A, &result);
  ck_assert_uint_eq(res, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sub_matrix) {
  matrix_t A, B, C, result, resultoriginal;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 3, &B);
  double numberA[] = {1,   2, 3, 4.8, 6,   8, 5.5, 7, -1, 2, 3,   4, -6,
                      8.5, 5, 7, -6,  8.4, 5, 7,   1, 0,  5, 1.2, 3};
  filling_arr(&A, numberA);
  double numberB[] = {1, 2, -3, 8, -5.5, 7, 3, -4, 6, 7, 6, -8.4};
  filling_arr(&B, numberB);
  // корректные матрицы одного размера
  int res = s21_sub_matrix(&A, &A, &result);
  ck_assert_int_eq(res, 0);
  s21_create_matrix(5, 5, &resultoriginal);
  double numberoriginalA[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  filling_arr(&resultoriginal, numberoriginalA);
  res = s21_eq_matrix(&result, &resultoriginal);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&resultoriginal);

  res = s21_sub_matrix(&B, &B, &result);
  ck_assert_int_eq(res, 0);
  s21_create_matrix(4, 3, &resultoriginal);
  double numberoriginalB[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  filling_arr(&resultoriginal, numberoriginalB);
  res = s21_eq_matrix(&result, &resultoriginal);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&resultoriginal);

  // корректные матрицы разных размеров
  res = s21_sub_matrix(&A, &B, &result);
  ck_assert_uint_eq(res, 2);
  s21_remove_matrix(&result);

  // некорректная матрица
  s21_create_matrix(0, -3, &C);
  res = s21_sub_matrix(&C, &C, &result);
  ck_assert_uint_eq(res, 1);
  res = s21_sub_matrix(&C, &A, &result);
  ck_assert_uint_eq(res, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_mult_number) {
  matrix_t A, B, C, result;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 3, &B);
  s21_create_matrix(5, 5, &C);
  double numberA[] = {1,   2, 3, 4.8, 6,   8, 5.5, 7, -1, 2, 3,   4, -6,
                      8.5, 5, 7, -6,  8.4, 5, 7,   1, 0,  5, 1.2, 3};
  double numberB[] = {1, 2, -3, 8, -5.5, 7, 3, -4, 6, 7, 6, -8.4};
  filling_arr(&A, numberA);
  filling_arr(&B, numberB);
  int res = s21_mult_number(&A, 0.00000001, &result);
  ck_assert_uint_eq(res, 0);
  res = s21_eq_matrix(&result, &C);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&C);

  res = s21_mult_number(&B, 3, &result);
  ck_assert_uint_eq(res, 0);
  double numberoriginB[] = {3, 6, -9, 24, -16.5, 21, 9, -12, 18, 21, 18, -25.2};
  filling_arr(&B, numberoriginB);
  res = s21_eq_matrix(&B, &result);
  ck_assert_uint_eq(res, 1);
  s21_remove_matrix(&result);

  s21_create_matrix(0, -3, &C);
  res = s21_mult_number(&C, 3, &result);
  ck_assert_uint_eq(res, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_matrix) {
  matrix_t A, B, C, D, result, original;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 3, &B);
  s21_create_matrix(3, 3, &D);
  double numberA[] = {1,   2, 3, 4.8, 6,   8, 5.5, 7, -1, 2, 3,   4, -6,
                      8.5, 5, 7, -6,  8.4, 5, 7,   1, 0,  5, 1.2, 3};
  filling_arr(&A, numberA);
  double numberB[] = {1, 2, -3, 8, -5.5, 7, 3, -4, 6, 7, 6, -8.4};
  filling_arr(&B, numberB);
  double numberD[] = {1, 2, -3, 8, -5.5, 7, 3, -4, 6};
  filling_arr(&D, numberD);
  // корректные матрицы одного размера
  int res = s21_mult_matrix(&A, &A, &result);
  ck_assert_int_eq(res, 0);
  s21_create_matrix(5, 5, &original);
  double numberoriginalA[] = {65.6, -3.8,  69.32, 59.5, 76.6,  68,    80.25,
                              22.1, 89.8,  93,    81.5, -47,   169.4, 7.9,
                              70.5, 26.2,  -15.4, 5.6,  144.4, 128,   27.4,
                              14.8, -1.92, 56.9,  48.4};
  filling_arr(&original, numberoriginalA);
  res = s21_eq_matrix(&result, &original);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&original);

  res = s21_mult_matrix(&B, &D, &result);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&D);
  s21_create_matrix(4, 3, &D);
  double numberoriginalD[] = {8,   3, -7, -15,  18.25, -20.5,
                              -11, 4, -1, 29.8, 14.6,  -29.4};
  filling_arr(&D, numberoriginalD);
  res = s21_eq_matrix(&D, &result);
  ck_assert_int_eq(res, 1);

  res = s21_mult_matrix(&B, &B, &result);
  ck_assert_int_eq(res, 2);

  // корректные матрицы разных размеров
  res = s21_mult_matrix(&A, &B, &result);
  ck_assert_uint_eq(res, 2);
  s21_remove_matrix(&result);

  // некорректная матрица
  s21_create_matrix(0, -3, &C);
  res = s21_mult_matrix(&C, &C, &result);
  ck_assert_uint_eq(res, 1);
  res = s21_mult_matrix(&C, &A, &result);
  ck_assert_uint_eq(res, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&D);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_transpose) {
  matrix_t A, B, C, result, resultorigin;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 3, &B);
  s21_create_matrix(5, -5, &C);
  double numberA[] = {1,   2, 3, 4.8, 6,   8, 5.5, 7, -1, 2, 3,   4, -6,
                      8.5, 5, 7, -6,  8.4, 5, 7,   1, 0,  5, 1.2, 3};
  double numberB[] = {1, 2, -3, 8, -5.5, 7, 3, -4, 6, 7, 6, -8.4};
  filling_arr(&A, numberA);
  filling_arr(&B, numberB);
  int res = s21_transpose(&A, &result);
  ck_assert_int_eq(res, 0);
  double numberoriginA[] = {1,   8, 3,   7,  1,   2, 5.5, 4, -6, 0, 3, 7, -6,
                            8.4, 5, 4.8, -1, 8.5, 5, 1.2, 6, 2,  5, 7, 3};
  s21_create_matrix(5, 5, &resultorigin);
  filling_arr(&resultorigin, numberoriginA);
  res = s21_eq_matrix(&result, &resultorigin);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&resultorigin);

  res = s21_transpose(&B, &result);
  ck_assert_int_eq(res, 0);
  double numberoriginB[] = {1, 8, 3, 7, 2, -5.5, -4, 6, -3, 7, 6, -8.4};
  s21_create_matrix(3, 4, &resultorigin);
  filling_arr(&resultorigin, numberoriginB);
  res = s21_eq_matrix(&result, &resultorigin);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&resultorigin);

  res = s21_transpose(&C, &result);
  ck_assert_int_eq(res, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_determinant) {
  matrix_t A, B, C, D, E, F;
  double result = 0, resultorigin = -2436.31;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 3, &B);
  s21_create_matrix(5, -5, &C);
  s21_create_matrix(5, 5, &D);
  s21_create_matrix(3, 3, &E);
  s21_create_matrix(1, 1, &F);
  double numberA[] = {1,   2, 3, 4.8, 6,   8, 5.5, 7, -1, 2, 3,   4, -6,
                      8.5, 5, 7, -6,  8.4, 5, 7,   1, 0,  5, 1.2, 3};
  double numberB[] = {1, 2, -3, 8, -5.5, 7, 3, -4, 6, 7, 6, -8.4};
  double numberE[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  filling_arr(&A, numberA);
  filling_arr(&B, numberB);
  filling_diag(&D);
  filling_arr(&E, numberE);
  filling_matrix(&F, 0.0000008);
  // корректная матрица
  int res = s21_determinant(&A, &result);
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq_tol(result, resultorigin, 0.0000001);

  // корректная матрица по диаг
  res = s21_determinant(&D, &result);
  resultorigin = 1;
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq_tol(result, resultorigin, 0.0000001);
  //  корректная матрица цифры по порядку
  res = s21_determinant(&E, &result);
  resultorigin = 0;
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq_tol(result, resultorigin, 0.0000001);
  // одиночная матрица
  res = s21_determinant(&F, &result);
  resultorigin = 0.0000008;
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq(result, resultorigin);
  // прямоугольная матрица
  res = s21_determinant(&B, &result);
  ck_assert_int_eq(res, 2);
  // некорректная матрица
  res = s21_determinant(&C, &result);
  ck_assert_int_eq(res, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&D);
  s21_remove_matrix(&E);
  s21_remove_matrix(&F);
}
END_TEST

START_TEST(test_s21_calc_complements) {
  matrix_t A, B, C, D, E, result, resultorigin;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 3, &B);
  s21_create_matrix(5, -5, &C);
  s21_create_matrix(3, 3, &D);
  s21_create_matrix(1, 1, &E);
  double numberA[] = {1,   2, 3, 4.8, 6,   8, 5.5, 7, -1, 2, 3,   4, -6,
                      8.5, 5, 7, -6,  8.4, 5, 7,   1, 0,  5, 1.2, 3};
  double numberB[] = {1, 2, -3, 8, -5.5, 7, 3, -4, 6, 7, 6, -8.4};
  double numberE[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  filling_arr(&A, numberA);
  filling_arr(&B, numberB);
  filling_matrix(&D, 5);
  filling_matrix(&E, 0.088);
  // корректная матрица
  int res = s21_calc_complements(&A, &result);
  ck_assert_int_eq(res, 0);
  s21_create_matrix(5, 5, &resultorigin);
  double numberoriginalA[] = {-160.79, 262.88,  1647.7,   3504.6,  -4094.41,
                              -228.72, -94.26,  171,      482,     -401.56,
                              69.46,   -212.96, -652.7,   -1697.3, 1743.6,
                              -216.95, 265.3,   270.85,   478.5,   -570.5,
                              864.51,  -727.02, -2953.55, -5618.2, 6069.59};
  filling_arr(&resultorigin, numberoriginalA);
  res = s21_eq_matrix(&result, &resultorigin);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&resultorigin);
  // корректная матрица из одинаковыых чисел
  res = s21_calc_complements(&D, &result);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&result);
  // одиночная матрица
  res = s21_calc_complements(&E, &result);
  ck_assert_int_eq(res, 2);
  s21_remove_matrix(&result);
  // прямоугольная матрица
  res = s21_calc_complements(&B, &result);
  ck_assert_int_eq(res, 2);
  // некоректная матрица
  res = s21_calc_complements(&C, &result);
  ck_assert_int_eq(res, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&D);
  s21_remove_matrix(&E);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_inverse_matrix) {
  matrix_t A, B, C, D, E, F, result, result1;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(4, 3, &B);
  s21_create_matrix(5, -5, &C);
  s21_create_matrix(5, 5, &D);
  s21_create_matrix(3, 3, &E);
  s21_create_matrix(1, 1, &F);
  double numberA[] = {1,   2, 3, 4.8, 6,   8, 5.5, 7, -1, 2, 3,   4, -6,
                      8.5, 5, 7, -6,  8.4, 5, 7,   1, 0,  5, 1.2, 3};
  double numberE[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  filling_arr(&A, numberA);
  filling_arr(&E, numberE);
  filling_diag(&D);
  filling_matrix(&F, 0.088);
  // корректная матрица det!=0
  int res = s21_inverse_matrix(&A, &result);
  ck_assert_int_eq(res, 0);
  s21_mult_matrix(&A, &result, &result1);
  s21_remove_matrix(&A);
  res = s21_eq_matrix(&D, &result1);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&D);
  s21_remove_matrix(&result1);

  // корректная матрица det==0
  res = s21_inverse_matrix(&E, &result);
  ck_assert_int_eq(res, 2);
  s21_remove_matrix(&result);

  // одиночная матрица
  res = s21_inverse_matrix(&F, &result);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&result);

  // прямоугольная
  res = s21_inverse_matrix(&B, &result);
  ck_assert_int_eq(res, 2);
  s21_remove_matrix(&result);

  // не корректная матрица
  res = s21_inverse_matrix(&C, &result);
  ck_assert_int_eq(res, 1);

  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&D);
  s21_remove_matrix(&E);
  s21_remove_matrix(&F);
  s21_remove_matrix(&result);
}
END_TEST

Suite *s21_matrix_suite(void) {
  Suite *suite;
  TCase *core;
  suite = suite_create("S21_matrix");
  core = tcase_create("Core");
  tcase_add_test(core, test_s21_create_matrix);
  tcase_add_test(core, test_s21_eq_matrix);
  tcase_add_test(core, test_s21_sum_matrix);
  tcase_add_test(core, test_s21_sub_matrix);
  tcase_add_test(core, test_s21_mult_number);
  tcase_add_test(core, test_s21_mult_matrix);
  tcase_add_test(core, test_s21_sub_matrix);
  tcase_add_test(core, test_s21_transpose);
  tcase_add_test(core, test_s21_determinant);
  tcase_add_test(core, test_s21_calc_complements);
  tcase_add_test(core, test_s21_inverse_matrix);

  suite_add_tcase(suite, core);

  return (suite);
}

int main(void) {
  int failed = 0;
  Suite *suite;

  SRunner *runner;

  suite = s21_matrix_suite();
  runner = srunner_create(suite);

  srunner_run_all(runner, CK_NORMAL);
  failed = srunner_ntests_failed(runner);
  srunner_free(runner);

  return (failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}

void filling_matrix(matrix_t *A, double number) {
  for (int i = 0; i < A->rows; i++) {
    for (int y = 0; y < A->columns; y++) {
      A->matrix[i][y] = number;
    }
  }
}

void filling_diag(matrix_t *A) {
  for (int i = 0, y = 0; i < A->rows; i++, y++) {
    A->matrix[i][y] = 1;
  }
}

void filling_arr(matrix_t *A, double *number) {
  int z = 0;
  for (int i = 0; i < A->rows; i++) {
    for (int y = 0; y < A->columns; y++) {
      A->matrix[i][y] = number[z];
      z++;
    }
  }
}
