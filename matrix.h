#pragma once
#include <vector>
const double TOL = 0.0000001;



template<typename T>
class Matrix
{
	// Custom Matrix class. Much of the code is taken directly from Antoine Savine's ptrMatrix class
	// Has an inversion method implemented using gaussian elimination (not the best method for inversion)
public:
	size_t ncols, nrows;
	std::vector<T> internal_vector;

	Matrix() : ncols(0), nrows(0) {};
	Matrix(const size_t rows, const size_t cols)
		: ncols(cols), nrows(rows), internal_vector(std::vector<T>(cols*rows)) {};
	Matrix(const size_t rows, const size_t cols, const std::vector<T>& vec)
		: ncols(cols), nrows(rows), internal_vector(vec) {};

	~Matrix() {};

		T* operator[] (const size_t row) { return &internal_vector[row*ncols]; }
	const T* operator[] (const size_t row) const { return &internal_vector[row*ncols]; }

	void resize(const size_t rows, const size_t cols)
	{
		ncols = cols;
		nrows = rows;
		internal_vector.resize(ncols*nrows);
	};

	void clear()
	{
		ncols = 0;
		nrows = 0;
		internal_vector.clear();
	}

	bool square() const { return ncols == nrows ? true : false; }

	void print() {
		std::cout << "Printing matrix with (rows x cols) " << nrows << " x " << ncols << std::endl;
		for (size_t m = 0; m < nrows; m++)
		{
			for (size_t n = 0; n < ncols; n++)
			{
				std::cout << internal_vector[m*ncols + n] << "  ";
			}
			std::cout << std::endl;
		}
	};

	void swap(Matrix<T>& rhs)
	{
		std::swap(internal_vector, rhs.internal_vector);
		std::swap(nrows, rhs.nrows);
		std::swap(ncols, rhs.ncols);
	}

	// Assign
	Matrix<T>& operator=(Matrix<T>& rhs)
	{
		Matrix<T> temp(rhs);
		swap(temp);
		return *this;
	}

	// Scalar operators
	Matrix<T>& operator+=(T& a);
	Matrix<T>& operator-=(T& a);
	Matrix<T>& operator*=(T a);

	Matrix<T> operator-() const;
	bool operator==(const Matrix& rhs) const;


	// Matrix multiplication
	Matrix<T> operator*(Matrix<T>& b);

	void transpose(); //also implemented is transpose(mat) which makes Matrix
	void inverse();   //also implement is function inverse(mat) which returns newmat
};

template<typename T>
inline bool Matrix<T>::operator==(const Matrix & rhs) const
{
	if (nrows != rhs.nrows || ncols != rhs.ncols) return false;
	for (size_t i = 0; i < nrows; ++i)
	{
		const double* ai = operator[](i);
		const double* bi = rhs[i];
		for (size_t j = 0; j < ncols; ++j)
		{
			if (fabs(ai[j] - bi[j]) > 1.0e-12) return false;
		}
	}
	return true;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator-() const
{
	Matrix<T> tmp(nrows, ncols);
	for (size_t i = 0; i < ncols*nrows; i++)
	{
		tmp.internal_vector[i] = -internal_vector[i];
	}
	return tmp;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator+=(T& a)
{
	for (size_t i = 0; i < ncols*nrows; i++)
	{
		internal_vector[i] = internal_vector[i] + a;
	}
	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator-=(T& a)
{
	for (size_t i = 0; i < ncols*nrows; i++)
	{
		internal_vector[i] = internal_vector[i] - a;
	}
	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator*=(T a)
{
	for (size_t i = 0; i < ncols*nrows; i++)
	{
		internal_vector[i] = internal_vector[i] * a;
	}
	return *this;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator*(Matrix<T>& b)
{
	if (ncols != b.nrows) {
		std::cout << "DIMENSIONALTY PROBLEMS MATRIX MULT stop alt" << std::endl;
	}

	const size_t rows = nrows, cols = b.ncols, n = ncols;

	Matrix<T> c(rows, cols);
	//  outermost loop on result rows
	for (size_t i = 0; i < rows; ++i)
	{
		//  loop on result columns
		for (size_t j = 0; j < cols; ++j)
		{
			//  compute dot product
			T res = 0.0;
			for (size_t k = 0; k < n; ++k)
			{
				res = res + internal_vector[i*ncols + k] * b[k][j];
			}
			c[i][j] = res;
		}   //  columns
	}   //  rows
	return c;
}

template<typename T>
inline void Matrix<T>::transpose()
{
	std::vector<T> tmp(ncols*nrows);
	for (size_t i = 0; i < ncols; i++)
	{
		for (size_t j = 0; j < nrows; j++)
		{
			tmp[nrows*i + j] = internal_vector[ncols*j + i];
		}
	}
	std::swap(ncols, nrows);
	internal_vector = tmp;
}

template<typename T>
Matrix<T> transpose(const Matrix<T>& mat) {
	Matrix<T> tmp = mat;
	tmp.transpose();
	return tmp;
}

template<typename T>
Matrix<T> identity(size_t n) {
	Matrix<T> tmp(n, n);

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (i == j) tmp[i][i] = 1.0;
			else tmp[i][j] = 0.0;
		}
	}

	return tmp;
}


template<typename T>
inline void Matrix<T>::inverse()
{
	Matrix<T> a(nrows, ncols);
	a.internal_vector = internal_vector;


	size_t n = ncols;
	size_t m = n;
	if (n != m)
		std::cout << "ERROR in .inverse(): not square matrix" << std::endl;

	size_t icol, irow;
	T big, dum, pivinv;

	Matrix<T> b = identity<T>(n);

	std::vector<size_t> indxc(n);
	std::vector<size_t> indxr(n);
	std::vector<size_t> ipiv(n);

	// Iterators
	size_t i, j, k, l, ll;

	for (j = 0; j < n; j++)
	{
		ipiv[j] = 0;
	}

	for (i = 0; i < n; i++)
	{
		big = 0.0;
		for (j = 0; j < n; j++)
		{
			if (ipiv[j] != 1)
				for (k = 0; k < n; k++)
				{
					if (ipiv[k] == 0) {
						if (abs(a[j][k]) >= big) {
							big = abs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}
		}
		ipiv[icol] = ipiv[icol] + 1;
		if (irow != icol)
		{
			for (l = 0; l < n; l++) { std::swap(a[irow][l], a[icol][l]); }
			for (l = 0; l < m; l++) { std::swap(b[irow][l], b[icol][l]); }
		}
		indxr[i] = irow;
		indxc[i] = icol;


		if (abs(a[icol][icol]) < TOL) {
			std::cout << "SINGULAR" << std::endl;
		}
		pivinv = 1.0 / a[icol][icol];
		for (l = 0; l < n; l++) { a[icol][l] = a[icol][l] * pivinv; }
		for (l = 0; l < m; l++) { b[icol][l] = b[icol][l] * pivinv; }
		for (ll = 0; ll < n; ll++)
		{
			if (ll != icol)
			{
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 0; l < n; l++) { a[ll][l] = a[ll][l] - a[icol][l] * dum; }
				for (l = 0; l < m; l++) { b[ll][l] = b[ll][l] - b[icol][l] * dum; }
			}
		}
	}
	swap(b);
}

template<typename T>
Matrix<T> inverse(const Matrix<T>& mat) {
	Matrix<T> tmp = mat;
	tmp.inverse();
	return tmp;
}

template<typename T>
Matrix<T> diag(const std::vector<T>& vec, const size_t rows, const size_t cols) {
	Matrix<T> tmp(rows, cols);
	for (size_t i = 0; i < rows*cols; i++) tmp.internal_vector[i] = 0;

	for (size_t i = 0; i < vec.size(); i++)
	{
		tmp[i][i] = vec[i];
	}
	return tmp;
}

template<typename T>
void diaginverse(Matrix<T>& m) {
	double TOL = 0.00001;
	for (size_t i = 0; i < m.ncols; i++)
	{
		if (abs(m[i][i]) > TOL) {
			m[i][i] = 1 / m[i][i];
		}
		else
		{
			m[i][i] = 0.0;
		}
	}
}

template<typename T>
void matrixwriter(std::string location, Matrix<T>& matrix) {
	std::ofstream myFile;
	myFile.open(location);

	for (size_t path = 0; path < matrix.nrows; path++)
	{
		for (size_t date = 0; date < matrix.ncols; date++)
		{
			myFile << value(matrix[path][date]) << ",";
		}
		myFile << "\n";
	}

	myFile.flush();
	myFile.close();

}

template<typename T>
Matrix<double> as_double(Matrix<T>& mat) {
	Matrix<double> out(mat.nrows, mat.ncols);

	out.internal_vector.resize(mat.internal_vector.size());
	for (size_t i = 0; i < mat.internal_vector.size(); i++)
	{
		out[i] = value(mat.internal_vector[i]);
	}
	return out;
}

template<typename T>
std::vector<double> as_double(std::vector<T>& vec) {
	std::vector<double> out(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		out[i] = value(vec[i]);
	}
	return out;
}