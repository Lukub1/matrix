CC=gcc -std=c11 
CFLAGS=-c -Wall -Wextra -Werror 
DEBUG = -g

all: clean s21_matrix.a test gcov_report

s21_matrix: s21_matrix.o
	@$(CC) $(DEBUG) s21_matrix.o -lm -o s21_matrix
	@./s21_matrix

s21_matrix.o: s21_matrix.c s21_matrix.h
	@$(CC) $(CFLAGS)  s21_matrix.c s21_matrix.h

s21_matrix.a: s21_matrix.o
	@ar -rcs s21_matrix.a s21_matrix.o
	@ranlib $@
	@cp $@ lib$@

clean:
	@rm -rf *.out *.o *.gch  log.txt logTest.txt s21_matrix
	@rm -rf coverage *.a *.gcda *gcno *.info Test *.cfg .clang-format

OS = $(shell uname)
ifeq ($(OS), Linux)
	LINUX_LIBS=-lsubunit -lrt -lpthread -lm
	CC += -D OS_LINUX -g
	CHECK_LEAKS=valgrind --leak-check=full -s --log-file=log.txt
	OPEN_GCOV=firefox coverage/src/index.html
	OUT_LEAKS= grep log.txt
else
	LINUX_LIBS=
	CC += -D OS_MAC
	CHECK_LEAKS=CK_FORK=no leaks --atExit --
	OPEN_GCOV=open coverage/index.html
endif

test:	clean style
	@$(CC) -g s21_matrix.c s21_matrix_test.c -lcheck $(LINUX_LIBS) -o ./Test 
	@$(CHECK_LEAKS) ./Test

style: 
	@cp ../materials/linters/.clang-format .clang-format
	@clang-format -n *.c *.h
	@rm -rf .clang-format

gcov_report: s21_matrix.c s21_matrix_test.c s21_matrix.h
	@$(CC) -fprofile-arcs -ftest-coverage s21_matrix_test.c s21_matrix.c -o Test -lcheck $(LINUX_LIBS) $(GCOV)
	@./Test
	@lcov -c -d . -o coverage.info
	@genhtml coverage.info -o coverage
	@$(OPEN_GCOV)

s21_matrix_test.o: s21_matrix_test.c s21_matrix.h
	@$(CC) $(CFLAGS) s21_matrix_test.c
