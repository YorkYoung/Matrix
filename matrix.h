#pragma once
#ifndef __OMEGA_MATRIX_H__
#define __OMEGA_MATRIX_H__

#if __cplusplus < 201103
#error This header file needs c++11 support.
#endif

#include <complex>
#include <valarray>
#include <exception>
#include <initializer_list>

namespace Omega
{
	using namespace std;

	class MatrixMissMatchedException : public exception
	{
	public:
		MatrixMissMatchedException(size_t rowSize1, size_t colSize1, size_t rowSize2, size_t colSize2) :
			RowSize1(rowSize1), ColSize1(colSize1),
			RowSize2(rowSize2), ColSize2(colSize2)
		{

		}

		const char* what() const noexcept override
		{
			return "Matrix miss matched!";
		}

		size_t RowSize1, ColSize1, RowSize2, ColSize2;
	};

	template <typename T>
	class Matrix
	{
	private:
		valarray<T> mElems;
		size_t mRowSize, mColSize;
	public:
		Matrix() : 
			mElems(), mRowSize(0), mColSize(0)
		{
			
		}

		Matrix(size_t rowSize, size_t colSize) : 
			mElems(static_cast<T>(0), rowSize*colSize),
			mRowSize(rowSize), mColSize(colSize)
		{

		}

		Matrix(initializer_list<initializer_list<T>> list) :
			mElems(), mRowSize(0), mColSize(0)
		{
			size_t rowSize = list.size();
			if (rowSize == 0) return;
			size_t colSize = list.begin()->size();
			this->mElems.resize(rowSize*colSize);
			this->mRowSize = rowSize;
			this->mColSize = colSize;

			size_t i = 0;
			for (const initializer_list<T>& rowList : list)
			{
				size_t j = 0;
				for (const T& elem : rowList)
				{
					(*this)[i][j] = elem;
					++j;
					if (j > colSize) break;
				}
				++i;
			}
		}

		Matrix(Matrix<T>&& mat) :
			mElems(move(mat.mElems)), 
			mRowSize(mat.mRowSize), mColSize(mat.mColSize)
		{

		}

		Matrix(const Matrix<T>& mat) :
			mElems(mat.mElems), 
			mRowSize(mat.mRowSize), mColSize(mat.mColSize)
		{

		}

		Matrix<T>& operator = (Matrix<T>&& mat)
		{
			this->mElems = move(mat.mElems);
			this->mRowSize = mat.mRowSize;
			this->mColSize = mat.mColSize;
			return *this;
		}

		Matrix<T>& operator = (const Matrix<T>& mat)
		{
			this->mElems = mat.mElems;
			this->mRowSize = mat.mRowSize;
			this->mColSize = mat.mColSize;
			return *this;
		}

		size_t Size() const
		{
			return this->mRowSize * this->mColSize;
		}

		size_t RowSize() const
		{
			return this->mRowSize;
		}

		size_t ColSize() const
		{
			return this->mColSize;
		}

		T* operator [] (int index)
		{
			return &(this->mElems[index*mColSize]);
		}

		const T* operator [] (int index) const
		{
			return &(this->mElems[index*mColSize]);
		}		

		Matrix<T> operator + () const
		{
			return *this;
		}

		Matrix<T> operator - () const
		{
			Matrix<T> mat2(this->mRowSize, this->mColSize);
			for (size_t i = 0; i < this->mRowSize; ++i)
			{
				for (size_t j = 0; j < this->mColSize; ++j)
				{
					mat2[i][j] = -(*this)[i][j];
				}
			}
			return move(mat2);
		}

		Matrix<T> operator * (T num) const
		{
			Matrix<T> mat2(this->mRowSize, this->mColSize);
			for (size_t i = 0; i < this->mRowSize; ++i)
			{
				for (size_t j = 0; j < this->mColSize; ++j)
				{
					mat2[i][j] = (*this)[i][j] * num;
				}
			}
			return move(mat2);
		}

		Matrix<T> operator / (T num) const
		{
			Matrix<T> mat2(this->mRowSize, this->mColSize);
			for (size_t i = 0; i < this->mRowSize; ++i)
			{
				for (size_t j = 0; j < this->mColSize; ++j)
				{
					mat2[i][j] = (*this)[i][j] / num;
				}
			}
			return move(mat2);
		}

		Matrix<T>& operator *= (T num) const
		{
			for (size_t i = 0; i < this->mRowSize; ++i)
			{
				for (size_t j = 0; j < this->mColSize; ++j)
				{
					(*this)[i][j] *= num;
				}
			}
			return *this;
		}

		Matrix<T>& operator /= (T num) const
		{
			for (size_t i = 0; i < this->mRowSize; ++i)
			{
				for (size_t j = 0; j < this->mColSize; ++j)
				{
					(*this)[i][j] /= num;
				}
			}
			return *this;
		}		

		Matrix<T>& operator += (const Matrix<T>& mat2)
		{
			if (this->mRowSize != mat2.mRowSize
				|| this->mColSize() != mat2.mColSize)
			{
				throw MatrixMissMatchedException(
					this->mRowSize, this->mColSize,
					mat2.mRowSize, mat2.mColSize
					);
			}

			for (size_t i = 0; i < mat2.mRowSize; ++i)
			{
				for (size_t j = 0; j < mat2.mColSize; ++j)
				{
					(*this)[i][j] += mat2[i][j];
				}
			}
			return *this;
		}

		Matrix<T>& operator -= (const Matrix<T>& mat2)
		{
			if (this->mRowSize != mat2.mRowSize
				|| this->mColSize() != mat2.mColSize)
			{
				throw MatrixMissMatchedException(
					this->mRowSize, this->mColSize,
					mat2.mRowSize, mat2.mColSize
					);
			}

			for (size_t i = 0; i < mat2.mRowSize; ++i)
			{
				for (size_t j = 0; j < mat2.mColSize; ++j)
				{
					(*this)[i][j] -= mat2[i][j];
				}
			}
			return *this;
		}

		Matrix<T> operator + (const Matrix<T>& mat2) const
		{
			if (this->mRowSize != mat2.mRowSize
				|| this->mColSize() != mat2.mColSize)
			{
				throw MatrixMissMatchedException(
					this->mRowSize, this->mColSize,
					mat2.mRowSize, mat2.mColSize
					);
			}

			Matrix<T> mat3(this->mRowSize, this->mColSize);
			for (size_t i = 0; i < this->mRowSize; ++i)
			{
				for (size_t j = 0; j < this->mColSize; ++j)
				{
					mat3[i][j] = (*this)[i][j] + mat2[i][j];
				}
			}
			return move(mat3);
		}

		Matrix<T> operator - (const Matrix<T>& mat2) const
		{
			if (this->mRowSize != mat2.mRowSize
				|| this->mColSize() != mat2.mColSize)
			{
				throw MatrixMissMatchedException(
					this->mRowSize, this->mColSize,
					mat2.mRowSize, mat2.mColSize
					);
			}

			Matrix<T> mat3(this->mRowSize, this->mColSize);
			for (size_t i = 0; i < this->mRowSize; ++i)
			{
				for (size_t j = 0; j < this->mColSize; ++j)
				{
					mat3[i][j] = (*this)[i][j] - mat2[i][j];
				}
			}
			return move(mat3);
		}

		Matrix<T> operator * (const Matrix<T>& mat2)
		{
			if (this->mColSize!= mat2.mRowSize)
			{
				throw MatrixMissMatchedException(
					this->mRowSize, this->mColSize,
					mat2.mRowSize, mat2.mColSize
					);
			}

			Matrix<T> mat3(this->mRowSize, mat2.mColSize);
			for (size_t i = 0; i < mat3.mRowSize; ++i)
			{
				for (size_t j = 0; j < mat3.mColSize; ++j)
				{
					mat3[i][j] = 0;
					for (size_t k = 0; k < this->mColSize; ++k)
					{
						mat3[i][j] += (*this)[i][k] * mat2[k][j];
					}
				}
			}
			return move(mat3);
		}

		Matrix<T>& operator *= (const Matrix<T>& mat2)
		{
			return (*this) = (*this) * mat2;
		}

		bool operator == (const Matrix<T>& mat2)
		{
			for (size_t i = 0; i < mat2.mRowSize; ++i)
			{
				for (size_t j = 0; j < mat2.mColSize; ++j)
				{
					if ((*this)[i][j] != mat2[i][j]) return false;
				}
			}
			return true;
		}

		bool operator != (const Matrix<T>& mat2)
		{
			for (size_t i = 0; i < mat2.mRowSize; ++i)
			{
				for (size_t j = 0; j < mat2.mColSize; ++j)
				{
					if ((*this)[i][j] != mat2[i][j]) return true;
				}
			}
			return false;
		}
	};

	template <typename T>
	Matrix<T> operator * (T num, const Matrix<T>& mat1)
	{
		Matrix<T> mat2(mat1.RowSize(), mat1.ColSize());
		for (size_t i = 0; i < mat1.RowSize(); ++i)
		{
			for (size_t j = 0; j < mat1.ColSize(); ++j)
			{
				mat2[i][j] = num * mat1[i][j];
			}
		}
		return move(mat2);
	}

	template <typename T>
	Matrix<T> trans(const Matrix<T>& mat1)
	{
		Matrix<T> mat2(mat1.ColSize(), mat1.RowSize());
		for (size_t i = 0; i < mat2.RowSize(); ++i)
		{
			for (size_t j = 0; j < mat2.ColSize(); ++j)
			{
				mat2[i][j] = mat1[j][i];
			}
		}
		return move(mat2);
	}

	template <typename T>
	Matrix<complex<T>> conj(const Matrix<complex<T>>& mat1)
	{
		Matrix<T> mat2(mat1.RowSize(), mat1.ColSize());
		for (size_t i = 0; i < mat2.RowSize(); ++i)
		{
			for (size_t j = 0; j < mat2.ColSize(); ++j)
			{
				mat2[i][j] = conj(mat1[i][j]);
			}
		}
		return move(mat2);
	}

	template <typename T>
	Matrix<complex<T>> HermiteConj(const Matrix<complex<T>>& mat1)
	{
		Matrix<T> mat2(mat1.ColSize(), mat1.RowSize());
		for (size_t i = 0; i < mat2.RowSize(); ++i)
		{
			for (size_t j = 0; j < mat2.ColSize(); ++j)
			{
				mat2[i][j] = conj(mat1[j][i]);
			}
		}
		return move(mat2);
	}

	class MatrixNotSquareException : public exception
	{
	public:
		MatrixNotSquareException(size_t rowSize, size_t colSize) :
			RowSize(rowSize), ColSize(colSize)
		{

		}

		const char* what() const noexcept override
		{
			return "Matrix is not square!";
		}

		size_t RowSize, ColSize;
	};

	class MatrixStrangeException : public exception
	{
	public:
		MatrixStrangeException(size_t rank, size_t rowSize, size_t colSize) :
			Rank(rank), RowSize(rowSize), ColSize(colSize)
		{

		}

		const char* what() const noexcept override
		{
			return "Strange matrix!";
		}

		size_t Rank, RowSize, ColSize;
	};

	template <typename T>
	Matrix<T> inv(const Matrix<T>& mat1)
	{
		typedef decltype(abs(T(0))) AbsT;
		if (mat1.RowSize() != mat1.ColSize())
		{
			throw MatrixNotSquareException(mat1.RowSize(), mat1.ColSize());
		}
		Matrix<T> mat2 = mat1.Clone();
		size_t size = mat2.RowSize();
		valarray<size_t> iMax(size), jMax(size);
		for (size_t k = 0; k < size; ++k)
		{
			AbsT maxAbs = abs(mat2[k][k]);
			iMax[k] = jMax[k] = k;
			for (size_t i = k; i < size; ++i)
			{
				for (size_t j = k; j < size; ++j)
				{
					AbsT tmpAbs = abs(mat2[i][j]);
					if (tmpAbs > maxAbs)
					{
						maxAbs = tmpAbs;
						iMax[k] = i;
						jMax[k] = j;
					}
				}
			}
			if (maxAbs == AbsT(0))
			{
				throw MatrixStrangeException(k, size, size);
			}
			if (iMax[k] != k) for (size_t j = 0; j < size; ++j)
			{
				//swap mat2[k][j] & mat2[iMax[k]][j]
				T tmp = mat2[iMax[k]][j];
				mat2[iMax[k]][j] = mat2[k][j];
				mat2[k][j] = tmp;

			}
			if (jMax[k] != k) for (size_t i = 0; i < size; ++i)
			{
				//swap mat2[i][k] & mat2[i][jMax[k]]
				T tmp = mat2[i][jMax[k]];
				mat2[i][jMax[k]] = mat2[i][k];
				mat2[i][k] = tmp;
			}
			mat2[k][k] = T(1.0) / mat2[k][k];
			for (size_t j = 0; j < size; ++j)
			{
				if (j != k) mat2[k][j] *= mat2[k][k];
			}
			for (size_t i = 0; i < size; ++i)
			{
				if (i != k) for (size_t j = 0; j < size; ++j)
				{
					if (j != k) mat2[i][j] -= mat2[i][k] * mat2[k][j];
				}
			}
			for (size_t i = 0; i < size; ++i)
			{
				if (i != k) mat2[i][k] = -mat2[i][k] * mat2[k][k];
			}
		}
		for (size_t l = 0, k = size - 1; l < size; ++l, --k)
		{
			if (jMax[k] != k) for (size_t j = 0; j < size; ++j)
			{
				//swap mat2[k][j] & mat2[jMax[k]][j]
				T tmp = mat2[jMax[k]][j];
				mat2[jMax[k]][j] = mat2[k][j];
				mat2[k][j] = tmp;
			}
			if (iMax[k] != k) for (size_t i = 0; i < size; ++i)
			{
				//swap mat2[i][k] & mat2[i][iMax[k]]
				T tmp = mat2[i][iMax[k]];
				mat2[i][iMax[k]] = mat2[i][k];
				mat2[i][k] = tmp;
			}
		}
		return move(mat2);
	}

	template <typename T>
	Matrix<T> ginv(const Matrix<T>& mat1)
	{
		size_t rowSize = mat1.ColSize(), colSize = mat1.RowSize();
		Matrix<T> mat2(rowSize, colSize);
		valarray<T> vectorB(colSize), vectorC(colSize), vectorD(rowSize);
		for (size_t k = 0; k < rowSize; ++k)
		{
			for (size_t i = 0; i < k; ++i)
			{
				vectorD[i] = 0;
				for (size_t j = 0; j < colSize; ++j)
				{
					vectorD[i] += mat2[i][j] * mat1[j][k];
				}
			}
			for (size_t i = 0; i < colSize; ++i)
			{
				vectorC[i] = mat1[i][k];
				for (size_t j = 0; j < k; ++j)
				{
					vectorC[i] -= mat1[i][j] * vectorD[j];
				}
			}
			T innerProduct = 0;
			T mod2 = 0;
			for (size_t i = 0; i < colSize; ++i)
			{
				innerProduct += vectorC[i] * vectorC[i];
				mod2 = abs(innerProduct);
			}
			if (mod2 == T(0))
			{
				innerProduct = 1;
				for (size_t i = 0; i < k; ++i)
				{
					innerProduct += vectorD[i] * vectorD[i];
				}
				for (size_t j = 0; j < colSize; ++j)
				{
					vectorB[j] = 0;
					for (size_t i = 0; i < k; ++i)
					{
						vectorB[j] += vectorD[i] * mat2[i][j];
					}
					vectorB[j] /= innerProduct;
				}
			}
			else
			{
				for (size_t j = 0; j < colSize; ++j)
				{
					vectorB[j] = vectorC[j] / innerProduct;
				}
			}
			for (size_t i = 0; i < k; ++i)
			{
				for (size_t j = 0; j < colSize; ++j)
				{
					mat2[i][j] -= vectorD[i] * vectorB[j];
				}
			}
			for (size_t j = 0; j < colSize; ++j)
			{
				mat2[k][j] = vectorB[j];
			}
		}
		return move(mat2);
	}

	template <typename T>
	Matrix<complex<T>> ginv(const Matrix<complex<T>>& mat1)
	{
		size_t rowSize = mat1.ColSize(), colSize = mat1.RowSize();		
		Matrix<complex<T>> mat2(rowSize, colSize);
		valarray<complex<T>> vectorB(colSize), vectorC(colSize), vectorD(rowSize);
		for (size_t k = 0; k < rowSize; ++k)
		{	
			for (size_t i = 0; i < k; ++i)
			{
				vectorD[i] = 0;
				for (size_t j = 0; j < colSize; ++j)
				{
					vectorD[i] += mat2[i][j] * mat1[j][k];
				}
			}
			for (size_t i = 0; i < colSize; ++i)
			{
				vectorC[i] = mat1[i][k];
				for (size_t  j = 0; j < k; ++j)
				{
					vectorC[i] -= mat1[i][j] * vectorD[j];
				}
			}
			complex<T> innerProduct(0, 0);
			T mod2 = 0;
			for (size_t i = 0; i < colSize; ++i)
			{
				innerProduct += conj(vectorC[i]) * vectorC[i];
				mod2 = abs(innerProduct);
			}
			if (mod2 == T(0))
			{
				innerProduct = complex<T>(1, 0);
				for (size_t i = 0; i < k; ++i)
				{
					innerProduct += conj(vectorD[i]) * vectorD[i];
				}
				for (size_t j = 0; j < colSize; ++j)
				{
					vectorB[j] = 0;
					for (size_t i = 0; i < k; ++i)
					{
						vectorB[j] += conj(vectorD[i]) * mat2[i][j];
					}
					vectorB[j] /= innerProduct;
				}
			}
			else
			{
				for (size_t j = 0; j < colSize; ++j)
				{
					vectorB[j] = conj(vectorC[j]) / innerProduct;
				}
			}
			for (size_t i = 0; i < k; ++i)
			{
				for (size_t j = 0; j < colSize; ++j)
				{
					mat2[i][j] -= vectorD[i] * vectorB[j];
				}
			}
			for (size_t j = 0; j < colSize; ++j)
			{
				mat2[k][j] = vectorB[j];
			}
		}
		return move(mat2);
	}

	template <typename T>
	Matrix<T> Triangulate(const Matrix<T>& mat1)
	{
		typedef decltype(abs(T(0))) AbsT;
		Matrix<T> mat2 = mat1;
		size_t rowSize = mat2.RowSize(), colSize = mat2.ColSize();
		size_t minSize = rowSize < colSize ? rowSize : colSize;		
		T sign = 1;
		for (size_t k = 0; k < minSize; ++k)
		{
			size_t iMax, jMax;
			AbsT maxAbs = 0, tmpAbs = 0;
			for (size_t i = k; i < rowSize; ++i) 
			{
				for (size_t j = k; j < colSize; ++j)
				{
					AbsT tmpAbs = abs(mat2[i][j]);
					if (tmpAbs > maxAbs)
					{
						maxAbs = tmpAbs;
						iMax = i;
						jMax = j;
					}
				}
			}
			if (maxAbs == AbsT(0)) break;
			if (iMax != k)
			{
				for (size_t j = 0; j < colSize; ++j)
				{
					//swap(mat2[k][j], mat2[iMax][j]);
					T tmp = mat2[k][j];
					mat2[k][j] = mat2[iMax][j];
					mat2[iMax][j] = tmp;
				}
				sign = -sign;
			}
			if (jMax != k)
			{
				for (size_t i = 0; i < rowSize; ++i)
				{
					//swap(mat2[i][k], mat2[i][jMax]);
					T tmp = mat2[i][k];
					mat2[i][k] = mat2[i][jMax];
					mat2[i][jMax] = tmp;
				}
				sign = -sign;
			}

			for (size_t i = k + 1; i < rowSize; ++i)
			{
				T quotient = mat2[i][k] / mat2[k][k];
				mat2[i][k] = 0;
				for (size_t j = k + 1; j < colSize; ++j)
				{
					mat2[i][j] -= quotient * mat2[k][j];
				}
			}
		}
		for (size_t i = minSize; i < rowSize; ++i)
		{
			for (size_t j = minSize; j < colSize; ++j)
			{
				mat2[i][j] = 0;
			}
		}
		mat2[0][0] *= sign;
		return move(mat2);
	}

	template <typename T>
	T det(const Matrix<T>& mat1)
	{
		size_t rowSize = mat1.RowSize(), colSize = mat1.ColSize();
		if (rowSize != colSize)
		{
			throw MatrixNotSquareException(rowSize, colSize);
		}
		Matrix<T> triangleMatrix = Triangulate(mat1);
		T determinant = 1;
		for (size_t i = 0; i < rowSize; ++i)
		{
			determinant *= triangleMatrix[i][i];
		}
		return determinant;
	}

	template <typename T>
	size_t rank(const Matrix<T>& mat1)
	{
		Matrix<T> triangleMatrix = Triangulate(mat1);
		size_t matRank = 0;
		size_t rowSize = mat1.RowSize(), colSize = mat1.ColSize();
		size_t minSize = rowSize < colSize ? rowSize : colSize;
		while (matRank < minSize && triangleMatrix[matRank][matRank] != 0)
		{
			++matRank;
		}
		return matRank;
	}
}

#endif // __OMEGA_MATRIX_H__
