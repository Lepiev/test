#include <iostream>
#include <fstream>

// Функция для освобождения памяти матрицы
void freeMatrix(double** matrix, int size) {
	for (int i = 0; i < size; ++i) {
		delete[] matrix[i];
	}
	delete[] matrix;
}

// Функция для умножения двух матриц
double** um(double** a, double** b, int size) {
	double** c = new double* [size];
	for (int i = 0; i < size; i++)
	{
		c[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			c[i][j] = 0;
			for (int k = 0; k < size; k++) {
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return c;
}

//заранее возводим матрицу в макс нужную степень и сохраняем её
void precomputePowers(double** matrix, int size, double*** powers, int maxExp) {
	powers[0] = new double* [size];
	for (int i = 0; i < size; ++i) {
		powers[0][i] = new double[size];
		for (int j = 0; j < size; ++j) {
			powers[0][i][j] = (i == j) ? 1.0 : 0.0; // Единичная матрица для степени 0
		}
	}

	powers[1] = new double* [size];
	for (int i = 0; i < size; ++i) {
		powers[1][i] = new double[size];
		for (int j = 0; j < size; ++j) {
			powers[1][i][j] = matrix[i][j]; // Исходная матрица для степени 1
		}
	}

	for (int exp = 2; exp <= maxExp; ++exp) {
		powers[exp] = um(powers[exp - 1], matrix, size);
	}
}

//след матрицы в степени
double trace(double*** powers, int size, int n) {
	double sum = 0;
	for (int i = 0; i < size; ++i) {
		sum += powers[n][i][i];
	}
	return sum;
}

//коэффициенты характеристического полинома с индексом k
double coef(double*** powers, int N, int k, double* coefCache, bool* isComputed) {
	if (k == 0) {
		return 1.0;
	}

	if (isComputed[k]) {
		return coefCache[k];
	}

	double result = -trace(powers, N, k);
	for (int i = 1; i < k; ++i) {
		result -= coef(powers, N, i, coefCache, isComputed) * trace(powers, N, k - i);
	}
	result /= k;

	coefCache[k] = result;
	isComputed[k] = true;

	return result;
}

//определитель методом Гауса
double det(double** matrix, int s) {
	// Создание копии матрицы
	double** tempMatrix = new double* [s];
	for (int i = 0; i < s; ++i) {
		tempMatrix[i] = new double[s];
		for (int j = 0; j < s; ++j) {
			tempMatrix[i][j] = matrix[i][j];
		}
	}
	double d = 1;
	for (int i = 0; i < s; ++i) {
		if (tempMatrix[i][i] == 0) {
			int k;
			for (k = i + 1; k < s; ++k) {
				if (tempMatrix[k][i] != 0) {
					// Меняем строки местами
					for (int j = 0; j < s; ++j) {
						double temp = tempMatrix[i][j];
						tempMatrix[i][j] = tempMatrix[k][j];
						tempMatrix[k][j] = temp;
					}
					d *= -1;
					break;
				}
			}
			if (k == s) {
				freeMatrix(tempMatrix, s);
				// всё ноль
				return 0;
			}
		}
		// Приведение матрицы к верхнетреугольной форме
		for (int k = i + 1; k < s; ++k) {
			double factor = tempMatrix[k][i] / tempMatrix[i][i];
			for (int j = i; j < s; ++j) {
				tempMatrix[k][j] -= factor * tempMatrix[i][j];
			}
		}
		d *= tempMatrix[i][i];
	}
	freeMatrix(tempMatrix, s);

	return d;
}

//сравнение строк
bool srav(const char* str1, const char* str2) {
	int i = 0;
	while (str1[i] != '\0' && str2[i] != '\0') {
		if (str1[i] != str2[i]) {
			return false;
		}
		i++;
	}
	return str1[i] == str2[i];
}

//обратить матрицу
void inverseMatrix(double** matrix, int size, double** inverse) {
	// Создаем единичную матрицу
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			inverse[i][j] = (i == j) ? 1.0 : 0.0;
		}
	}

	// Применяем метод Гаусса-Жордана
	for (int i = 0; i < size; ++i) {
		// Обеспечиваем, чтобы элемент matrix[i][i] был ненулевым
		if (matrix[i][i] == 0.0) {
			// Ищем ненулевой элемент в текущем столбце
			for (int j = i + 1; j < size; ++j) {
				if (matrix[j][i] != 0.0) {
					// Меняем строки местами в обеих матрицах
					for (int k = 0; k < size; ++k) {
						std::swap(matrix[i][k], matrix[j][k]);
						std::swap(inverse[i][k], inverse[j][k]);
					}
					break;
				}
			}
		}

		// Нормализуем строку так, чтобы элемент matrix[i][i] был равен 1
		double factor = matrix[i][i];
		for (int j = 0; j < size; ++j) {
			matrix[i][j] /= factor;
			inverse[i][j] /= factor;
		}

		// Обнуляем все элементы в текущем столбце
		for (int j = 0; j < size; ++j) {
			if (j != i) {
				factor = matrix[j][i];
				for (int k = 0; k < size; ++k) {
					matrix[j][k] -= factor * matrix[i][k];
					inverse[j][k] -= factor * inverse[i][k];
				}
			}
		}
	}
}


int main(int argc, char* argv[]) {

	// Чтение аргументов
	const char* mode = argv[1];
	const char* inputFilename = argv[2];
	const char* outputFilename = argv[3];


	// Открываем файл для чтения матрицы
	std::fstream inputFile(inputFilename);
	if (!inputFile) {
		return 0;
	}

	// Чтение размера матрицы
	int rows;
	int cols;
	inputFile >> rows >> cols;


	// Выделение памяти под матрицу и чтение матрицы
	double** matrix = new double* [rows];
	for (int i = 0; i < rows; ++i) {
		matrix[i] = new double[cols];
		for (int j = 0; j < cols; ++j) {
			inputFile >> matrix[i][j];
		}
	}
	inputFile.close();

	// Открываем файл для записи результата
	std::ofstream outputFile(outputFilename);
	if (!outputFile) {
		freeMatrix(matrix, rows); // Освобождение памяти матрицы
		return 0;
	}

	// Проверка на квадратную матрицу
	if (rows != cols) {
		freeMatrix(matrix, rows); // Освобождение памяти матрицы
		outputFile << "error";
		outputFile.close();
		return 0;
	}

	// Выполнение расчета в зависимости от режима
	if (srav(mode, "-det")) {
		double determinant = det(matrix, rows);
		outputFile << determinant << "\n";
	}
	else if (srav(mode, "-poly")) {
		//кэш
		double* coefCache = new double[rows + 1];
		bool* isComputed = new bool[rows + 1] {false};


		// чтобы раз за разом не возводить матрицу в степень
		double*** powers = new double** [rows + 1];
		precomputePowers(matrix, rows, powers, rows);

		outputFile << 1 << " " << rows + 1 << std::endl;
		for (int i = 0; i <= rows; ++i) {
			double coefficient = coef(powers, rows, i, coefCache, isComputed);
			outputFile << coefficient << " ";
		}

		//удаляем макс степень
		for (int i = 0; i <= rows; ++i) {
			freeMatrix(powers[i], rows);
		}
		delete[] powers;
		delete[] coefCache;
		delete[] isComputed;
	}
	else if (srav(mode, "-inv")) {
		double determinant = det(matrix, rows);
		if (determinant == 0) {
			outputFile << "singular";
		}
		else {
			double** inv = new double* [rows];
			for (int i = 0; i < rows; ++i) {
				inv[i] = new double[rows];
			}
			inverseMatrix(matrix, rows, inv);
			outputFile << rows << " " << rows << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < rows; ++j) {
					outputFile << inv[i][j] << " ";
				}
				outputFile << std::endl;
			}
			freeMatrix(inv, rows);
		}
	}

	

	else {
		freeMatrix(matrix, rows); // Освобождение памяти матрицы
		outputFile.close();

		return 0;
	}

	// Очистка и закрытие файла
	freeMatrix(matrix, rows); // Освобождение памяти матрицы
	outputFile.close();


	return 0;
}
