# Функция для вычисления определителя 2x2 матрицы
def det2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


# Функция для получения минора матрицы
def minor(matrix, i, j):
    tmp = [row[:j] + row[j + 1:] for row in (matrix[:i] + matrix[i + 1:])]
    return tmp


# Рекурсивная функция для вычисления определителя NxN матрицы
def determinant(matrix):
    size = len(matrix)
    if size == 2:
        return det2(matrix)

    det = 0
    for j in range(size):
        det += ((-1) ** j) * matrix[0][j] * determinant(minor(matrix, 0, j))

    return det


# Функция для транспонирования матрицы
def transposeMatrix(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]


# Функция для получения обратной матрицы
def getMatrixInverse(matrix):
    det = determinant(matrix)
    if det == 0:
        raise ValueError("Матрица вырождена, обратной матрицы не существует")

    cofactors = []
    for i in range(len(matrix)):
        cofactorRow = []
        for j in range(len(matrix)):
            minorMtx = minor(matrix, i, j)
            cofactorRow.append(((-1) ** (i + j)) * determinant(minorMtx))
        cofactors.append(cofactorRow)

    adjugate = transposeMatrix(cofactors)
    for i in range(len(adjugate)):
        for j in range(len(adjugate)):
            adjugate[i][j] /= det

    return adjugate


# Функция для приведения матрицы к треугольному виду с выбором главного элемента
def make_triangle_pivot(matrix):
    for nrow in range(len(matrix)):
        pivot = max(range(nrow, len(matrix)), key=lambda i: abs(matrix[i][nrow]))
        matrix[nrow], matrix[pivot] = matrix[pivot], matrix[nrow]
        diag = matrix[nrow][nrow]
        if abs(diag) < 1e-10:
            raise ValueError("Матрица несовместна")
        matrix[nrow] = [x / diag for x in matrix[nrow]]
        for i in range(nrow + 1, len(matrix)):
            factor = matrix[i][nrow]
            matrix[i] = [x - factor * y for x, y in zip(matrix[i], matrix[nrow])]


# Функция для приведения матрицы к единичной
def make_identity(matrix):
    for nrow in range(len(matrix) - 1, -1, -1):
        for upper_row in matrix[:nrow]:
            factor = upper_row[nrow]
            upper_row[-1] -= factor * matrix[nrow][-1]
            upper_row[nrow] = 0


# Функция для решения системы линейных уравнений методом Гаусса
def gauss_solve(A, b=None):
    if b is not None:
        A = [row + [val] for row, val in zip(A, b)]
    make_triangle_pivot(A)
    make_identity(A)
    return [row[-1] for row in A]


# Исходная матрица A и вектор b
myMatrix = [[1, 2, -2, 6],
            [-3, -5, 14, 13],
            [1, 2, -2, -2],
            [-2, -4, 5, 10]]

myB = [24, 41, 0, 20]

# Решение системы уравнений и вычисление определителя и обратной матрицы
result = gauss_solve(myMatrix, myB)
det_A = determinant(myMatrix)
inverse_A = getMatrixInverse(myMatrix)

# Вывод результатов
print("Определитель матрицы A:", det_A)
print()
print("Обратная матрица A:")
for row in inverse_A:
    print(row)
print()
print("Решение системы уравнений:")
print(result)
