#pragma once
#ifndef _NUMERIC_MATRIX_
#define _NUMERIC_MATRIX_

#include "array.h"
#include "complex.h"

namespace numeric
{
	namespace algebra
	{
#pragma region _Matrix_
		template <typename value_type>
		class matrix : public array_base<value_type>
		{
		public:
			typedef array_base<value_type> my_base;
			typedef matrix<value_type> my_type;
			typedef value_type value_type;
			typedef size_t size_type;
			typedef ptrdiff_t difference_type;
			typedef value_type* pointer;
			typedef const value_type* const_pointer;
			typedef value_type& reference;
			typedef const value_type& const_reference;

			typedef array_iterator<value_type> iterator;
			typedef const_array_iterator<value_type> const_iterator;

			matrix() :
			my_base(),
			row_size(0),
			column_size(0)
			{ }

			matrix(size_type var_height, size_type var_width, value_type value = value_type()) :
			my_base(var_height * var_width, value),
			row_size(var_height),
			column_size(var_width)
			{ }

			matrix(const my_type& var_matrix) :
			my_base(var_matrix),
			row_size(var_matrix.height()),
			column_size(var_matrix.width())
			{ }

			template <typename another_type>
			matrix(const matrix<another_type>& var_matrix) :
			my_base(var_matrix), 
			row_size(var_matrix.height()), 
			column_size(var_matrix.width())
			{ }

			size_type size() const
			{
				return row_size * column_size;
			}

			void clear()
			{
				row_size = 0;
				column_size = 0;
			}

			size_type height() const
			{
				return row_size;
			}

			size_type width() const
			{
				return column_size;
			}

			void set_size(size_type var_height, size_type var_width)
			{
				size_type new_size = var_height * var_width;
				size_type move_size = std::min(row_size, var_height);
				if(new_size > capacity()) resize(new_size);				
				if(column_size > var_width)
				{
					for(size_t i = 1, j = 0; i < move_size; ++i) for(j = 0; j < var_width; ++j)
					{
						array_head[i * var_width + j] = array_head[i * column_size + j];
					}
				}
				if(column_size < var_width)
				{
					for(size_t i = move_size - 1, j = 0; i > 0; --i) for(j = 0; j < var_width; ++j)
					{
						array_head[i * var_width + j] = array_head[i * column_size + j];
					}
				}
				row_size = var_height;
				column_size = var_width;
			}

			bool empty() const
			{
				return row_size == 0 || column_size == 0; 
			}

			pointer operator [] (int index)
			{
				return array_head + index * column_size;
			}

			const_pointer operator [] (int index) const
			{
				return array_head + index * column_size;
			}

			my_type& operator = (const my_type& var_matrix)
			{
				copy_data(var_matrix);
				row_size = var_matrix.row_size;
				column_size = var_matrix.column_size;
				return *this;
			}

			template <typename another_type>
			my_type& operator = (const matrix<another_type>& var_matrix)
			{
				copy_data(var_matrix);
				row_size = var_matrix.row_size;
				column_size = var_matrix.column_size;
				return *this;
			}

			void swap(my_type& var_matrix)
			{
				my_base::swap(var_matrix);
				std::swap(row_size, var_matrix.row_size);
				std::swap(column_size, var_matrix.column_size);
			}

			void swap(my_type &&  var_matrix)
			{
				my_base::swap((my_base && )(var_array));
				std::swap(row_size, var_matrix.row_size);
				std::swap(column_size, var_matrix.column_size);
			}

			template <typename another_type>
			bool operator == (const matrix<another_type>& var_matirx) const
			{
				return equals(var_matirx);
			}

			template <typename another_type>
			bool operator != (const matrix<another_type>& var_matirx) const
			{
				return !equals(var_matirx);
			}

		protected:
			template <typename another_type>
			bool equals(const matrix<another_type>& var_matirx) const
			{
				if(row_size == var_matirx.height()  &&  column_size = var_matirx.width())
				{	
					const_iterator my_iter = begin();
					const_iterator var_iter = var_matirx.begin(), var_end = var_matirx.end();
					for(; *my_iter == *var_iter  &&  var_iter < var_end; ++my_iter, ++var_iter);
					return var_iter == var_end;
				}
				else return false;
			}

		private:
			size_type row_size, column_size;
		};

		template <typename value_type>
		inline void swap(matrix<value_type>& left, matrix<value_type>& right)
		{
			left.swap(right);
		}

		template <typename value_type>
		inline void swap(matrix<value_type> &&  left, matrix<value_type>& right)
		{
			right.swap((matrix<value_type> && )(left));
		}

		template <typename value_type>
		inline void swap(matrix<value_type>& left, matrix<value_type> &&  right)
		{
			left.swap((matrix<value_type> && )(right));
		}
#pragma endregion

#pragma region _Matrix_Operation_
		template <typename value_type>
		inline matrix<value_type> operator + (const matrix<value_type>& var_matrix)
		{
			return var_matrix;
		}

		template <typename value_type>
		matrix<value_type> operator - (const matrix<value_type>& var_matrix)
		{
			size_t height = var_matrix.height(), width = var_matrix.width(); 
			matrix<value_type> result(height, width);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
 			{
				result[i][j] = -var_matrix[i][j];
			}
			return result;
		}

		template <typename value_type>
		matrix<value_type> trans(const matrix<value_type>& var_matrix)
		{
			size_t height = var_matrix.height(), width = var_matrix.width(); 
			matrix<value_type> result(width, height);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				result[i][j] = -var_matrix[j][i];
			}
			return result;
		}

		template <typename value_type>
		matrix<value_type> conj(const matrix<value_type>& var_matrix)
		{
			size_t height = var_matrix.height(), width = var_matrix.width(); 
			matrix<value_type> result(height, width);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				result[i][j] = conj(var_matrix[i][j]);
			}
			return result;
		}

		template <typename value_type>
		matrix<value_type> tconj(const matrix<value_type>& var_matrix)
		{
			size_t height = var_matrix.height(), width = var_matrix.width(); 
			matrix<value_type> result(width, height);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				result[i][j] = conj(var_matrix[j][i]);
			}
			return result;
		}

		template <typename value_type, typename another_type>
		matrix<value_type> operator + (
			const matrix<value_type>& left_matrix,
			const matrix<another_type>& right_matrix)
		{
			size_t height = left_matrix.height(), width = left_matrix.width(); 
			if(height != right_matrix.height() || width != right_matrix.width())
			{
				throw exception("Matrix type mismatched!")
			}
			matrix<value_type> result(width, height);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				result[i][j] = left_matrix[i][j] + right_matrix[i][j];
			}
			return result;
		}

		template <typename value_type, typename another_type>
		matrix<value_type>& operator += (
			matrix<value_type>& A,
			const matrix<another_type>& B)
		{
			size_t height = left_matrix.height(), width = left_matrix.width(); 
			if(height != right_matrix.height() || width != right_matrix.width())
			{
				throw exception("Matrix type mismatched!")
			}
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				left_matrix[i][j] += right_matrix[i][j];
			}
			return left_matrix;
		}

		template <typename value_type, typename another_type>
		matrix<value_type> operator - (
			const matrix<value_type>& left_matrix,
			const matrix<another_type>& right_matrix)
		{
			size_t height = left_matrix.height(), width = left_matrix.width(); 
			if(height  != right_matrix.height() || width  != right_matrix.width())
			{
				throw exception("Matrix type mismatched!")
			}
			matrix<value_type> result(width, height);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				result[i][j] = left_matrix[i][j] - right_matrix[i][j];
			}
			return result;
		}

		template <typename value_type, typename another_type>
		matrix<value_type>& operator -= (
			matrix<value_type>& A,
			const matrix<another_type>& B)
		{
			size_t height = left_matrix.height(), width = left_matrix.width(); 
			if(height != right_matrix.height() || width  != right_matrix.width())
			{
				throw exception("Matrix type mismatched!")
			}
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				left_matrix[i][j] -= right_matrix[i][j];
			}
			return left_matrix;
		}

		template <typename value_type, typename another_type>
		matrix<value_type> operator * (
			const matrix<value_type>& left_matrix,
			const matrix<another_type>& right_matrix)
		{
			size_t height = left_matrix.height(), width = right_matrix.width();
			size_t multi_size = left_matrix.width();
			if(multi_size  != right_matrix.height())
			{
				throw exception("Matrix type mismatched!")
			}
			matrix<value_type> result(width, height);
			for(size_t i = 0, j = 0 , k = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				for(result[i][j] = 0, k = 0;  k < multi_size; ++k)
				{
					result[i][j] += left_matrix[i][k] * right_matrix[k][j];
				}
			}
			return result;
		}

		template <typename value_type, typename another_type>
		matrix<value_type>& operator *= (
			matrix<value_type>& left_matrix,
			const matrix<another_type>& right_matrix)
		{
			return left_matrix = left_matrix * right_matrix;
		}

		template <typename value_type, typename another_type>
		matrix<another_type> operator * (
			const value_type& var,
			const matrix<another_type>& var_matrix)
		{
			size_t height = var_matrix.height(), width = var_matrix.width(); 
			matrix<value_type> result(width, height);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				result[i][j] = var * var_matrix[i][j];
			}
			return result;
		}

		template <typename value_type, typename another_type>
		matrix<value_type> operator * (
			const matrix<value_type>& var_matrix,
			const another_type& var)
		{
			size_t height = var_matrix.height(), width = var_matrix.width(); 
			matrix<value_type> result(width, height);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				result[i][j] = var_matrix[i][j] * var;
			}
			return result;
		}

		template <typename value_type, typename another_type>
		matrix<value_type>& operator *= (
			const matrix<value_type>& var_matrix,
			const another_type& var)
		{
			size_t height = var_matrix.height(), width = var_matrix.width(); 
			matrix<value_type> result(width, height);
			for(size_t i = 0, j = 0;  i < height; ++i) for(j = 0; j < width; ++j)
			{
				var_matrix[i][j] *= var;
			}
			return var_matrix;
		}

		template <typename value_type>
		matrix<value_type> operator ^ (
			const matrix<value_type>& var_matrix, 
			int index)
		{
			size_t var_size = var_matrix.height();
			if(var_size != var_matrix.width()) throw ("The matrix is not square!");
			matrix<value_type> result(var_size, var_size), power;
			uint abs_index;
			for(size_t i = 0; i < var_size; ++i)
			{
				result[i][i] = 1;
			}			
			if(index >= 0)
			{
				power = var_matrix;
				abs_index = index;
			}
			else
			{
				power = inv(var_matrix);
				abs_index = -index;
			}
			for(;abs_index != 0; abs_index >>= 1)
			{
				if(abs_index & 1 == 1)
				{
					result *= power;
				}
				power *= power;
			}
			return result;
		}

		template <typename value_type, typename another_type>
		matrix<value_type> operator ^=  (
			matrix<value_type>& var_matrix, 
			int index)
		{
			size_t var_size = var_matrix.height();
			if(var_size != var_matrix.width()) throw ("The matrix is not square!");
			matrix<value_type> power;
			uint abs_index;
			for(size_t i = 0, j = 0; i < var_size; ++i) 
			{
				for(j = 0; j < i; ++j)
				{
					var_matrix[i][j] = 0;
				}
				var_matrix[i][j] = 1;
				for(++j; j < var_size; ++j)
				{
					var_matrix[i][j] = 0;
				}
			}			
			if(index >= 0)
			{
				power = var_matrix;
				abs_index = index;
			}
			else
			{
				power = inv(var_matrix);
				abs_index = -index;
			}
			for(;abs_index != 0; abs_index >>= 1)
			{
				if(abs_index & 1 == 1)
				{
					var_matrix *= power;
				}
				power *= power;
			}
			return var_matrix;
		}
#pragma endregion

#pragma region _Matrix_Inverse_
		template <typename value_type>
		matrix<value_type> inv(const matrix<value_type>& var_matrix)
		{
			matrix<value_type> result = var_matrix;
			size_t i, j, k;
			size_t var_size = result.height();
			ldouble max_abs, tmp_abs;
			array_list<size_t> i_max(var_size), j_max(var_size);
			if(var_size != result.width()) throw exception("The matrix is not square!");
			for(k = 0; k < var_size; ++k)
			{
				for(max_abs = 0, i = k; i < var_size; ++i) for(j = k; j < var_size; ++j)
				{
					tmp_abs = abs(result[i][j]);
					if(tmp_abs > max_abs)
					{
						max_abs = tmp_abs;
						i_max[k] = i;
						j_max[k] = j;
					}
				}
				if(max_abs + 1.0 == 1.0) throw exception("The matrix is singular!");
				if(i_max[k] != k) for(j = 0; j < var_size; ++j)
				{
					swap(result[k][j], result[i_max[k]][j]);
				}
				if(j_max[k] != k) for(i = 0; i < var_size; ++i)
				{
					swap(result[i][k], result[i][j_max[k]]);
				}
				result[k][k] = 1.0 / result[k][k];
				for(j = 0; j < var_size; ++j)
				{
					if(j != k) result[k][j] *= result[k][k];
				}
				for(i = 0; i < var_size; ++i)
				{
					if(i !=k) for(j=0; j < var_size; ++j)
					{  
						if(j != k) result[i][j] -= result[i][k] * result[k][j];
					}
				}
				for(i = 0; i < var_size; ++i) 
				{
					if(i != k) result[i][k] = -result[i][k] * result[k][k];
				}
			}
			for(k = var_size - 1; k != numeric::end; --k)
			{
				if(j_max[k] != k) for(j = 0; j < var_size; ++j)
				{ 
					swap(result[k][j], result[j_max[k]][j]);
				}
				if(i_max[k] != k) for(i = 0; i < var_size; ++i)
				{ 
					swap(result[i][k], result[i][i_max[k]]);
				}
			}
			return result;
		}

		template <typename value_type>
		matrix<value_type> ginv(const matrix<value_type>& var_matrix)
		{
			size_t i, j, k;
			size_t height = var_matrix.height(), width = var_matrix.width();
			value_type inner_product;
			ldouble mod2;
			matrix<value_type> result(width, height);
			array_list<value_type> vector_b(height), vector_c(height), vector_d(width);
			for(k = 0; k < width; ++k)
			{
				for(i = 0; i < k; ++i) for(vector_d[i] = 0, j = 0; j < height; ++j)
				{
					vector_d[i] += result[i][j] * var_matrix[j][k];
				}
				for(i = 0; i < height; ++i) for(vector_c[i] = var_matrix[i][k], j = 0; j < k; ++j) 
				{
					vector_c[i] -= var_matrix[i][j] * vector_d[j];
				}
				for(inner_product = 0, i = 0; i < height; ++i) 
				{
					inner_product += conj(vector_c[i]) * vector_c[i];
					mod2 = abs(inner_product);
				}
				if(mod2 + 1.0 == 1.0)
				{
					for(inner_product = 1, i = 0; i < k; ++i)
					{
						inner_product += conj(vector_d[i]) * vector_d[i];
					}
					for(j = 0; j < height; ++j)
					{
						for(vector_b[j]=0, i=0; i<k;++i) 
						{
							vector_b[j] += conj(vector_d[i]) * result[i][j];
						}
						vector_b[j] /= inner_product;
					}
				}
				else
				{
					for(j = 0; j < height; ++j)
					{
						vector_b[j] = conj(vector_c[j]) / inner_product;
					}
				}
				for(i = 0; i < k; ++i) for(j = 0; j < height; ++j) 
				{
					result[i][j] -= vector_d[i] * vector_b[j];
				}
				for(j = 0; j < height; ++j) 
				{
					result[k][j] = vector_b[j];
				}
			}
			return result;
		}
#pragma endregion

#pragma region _Matrix_Transform_
		template <typename value_type>
		matrix<value_type> triangulate(const matrix<value_type>& var_matrix)
		{
			matrix<value_type> result = var_matrix;
			size_t i, j, k, i_max, j_max;
			size_t height = result.height(), width = result.width();
			size_t min_size = std::min(height, width);
			ldouble max_abs, tmp_abs;
			value_type quotient;
			int sign = 1;
			for(k = 0; k < min_size; ++k)
			{
				for(max_abs = 0, i = k; i < height; ++i) for(j = k; j < width; ++j)
				{
					tmp_abs = abs(result[i][j]);
					if(tmp_abs > max_abs)
					{
						max_abs = tmp_abs;
						i_max = i;
						j_max = j;
					}
				}
				if(max_abs + 1.0 == 1.0) break;
				if(i_max != k)
				{
					for(j = 0; j < width; ++j)
					{
						swap(result[k][j], result[i_max][j]);
					}
					sign *= -1;
				}
				if(j_max != k)
				{
					for(i = 0; i < height; ++i)
					{
						swap(result[i][k], result[i][j_max]);
					}
					sign *= -1;
				}
				for(i = k + 1 ; i < height; ++i)
				{
					quotient = result[i][k] / result[k][k];
					result[i][k] = 0; 
					for(j = k + 1; j < width; ++j)
					{
						result[i][j] -= quotient * result[k][j];
					}
				}
			}
			for(i = k; i < height; ++i) for(j = k; j < width; ++j)
			{
				result[i][j] = 0; 
			}
			result[0][0] *= sign;
			return result;
		}

		template <typename value_type>
		matrix<value_type> normalize(const matrix<value_type>& var_matrix)
		{
			matrix<value_type> result = var_matrix;
			size_t i, j, k, l, i_max;
			size_t height = result.height(), width = result.width();
			ldouble max_abs, tmp_abs;
			value_type quotient;
			for(l = 0, k = 0; k < height && l < width; ++l)
			{
				for(max_abs = 0, i = k; i < height; ++i)
				{
					tmp_abs = abs(result[i][l]);
					if(tmp_abs > max_abs)
					{
						max_abs = tmp_abs;
						i_max = i;
					}
				}
				if(max_abs + 1.0 == 1.0) continue;
				if(i_max != k)
				{
					for(j = 0; j < width; ++j)
					{ 
						swap(result[k][j], result[i_max][j]);
					}
				}
				for(i = k + 1; i < height; ++i)
				{
					quotient = result[i][l] / result[k][l];
					result[i][l] = 0;
					for(j = l + 1; j < width; ++j)
						result[i][j] -= quotient * result[k][j];
				} 
				++k;
			}
			for(i = 0, j = 0; i < height; ++i)
			{
				for(; j < width && result[i][j] == 0;  ++j);
				if(j < width)
				{
					quotient = result[i][j];
					result[i][j] = 1;
					for(k = j + 1; k < width; ++k)
					{
						result[i][k] /= quotient;
					}
					for(l = 0; l < i; ++l)
					{
						quotient = result[l][j];
						result[l][j] = 0; 
						for(k = j + 1; k < width; ++k)
						{
							result[i_max][k]- = quotient * result[i][k];
						}
					}
				}
				else break;
			}
			return result;
		}

		template <typename value_type>
		matrix<value_type> linear_solve(const matrix<value_type>& var_matrix)
		{
			matrix<value_type> tmp_matrix = var_matrix, result;
			size_t i, j, k, i_max, j_max;
			size_t height = tmp_matrix.height(), width = tmp_matrix.width();
			size_t min_szie = std::min(height, width);
			ldouble max_abs, tmp_abs;
			value_type quotient;
			array_list<size_t> index_array(width);
			for(i = 0; i < width; ++i)
			{
				index_array[i] = i;
			}
			for(k = 0; k < min_szie; ++k)
			{
				for(max_abs = 0, i = k; i < height; ++i) for(j = k; j < width; ++j)
				{
					tmp_abs = abs(tmp_matrix[i][j]);
					if(tmp_abs > max_abs)
					{
						max_abs = tmp_abs;
						i_max = i;
						j_max = j;
					}
				}
				if(max_abs + 1.0 == 1.0) break;
				if(i_max != k)
				{
					for(j = 0; j < width; ++j)
					{
						swap(tmp_matrix[k][j], tmp_matrix[i_max][j]);
					}
				}
				if(j_max != k)
				{
					for(i = 0; i < height; ++i)
					{
						swap(tmp_matrix[i][k], tmp_matrix[i][j_max]);
					}
					swap(index_array[k], index_array[j_max]);
				}
				for(i = k + 1; i < height; ++i)
				{
					quotient = tmp_matrix[i][k] / tmp_matrix[k][k];
					tmp_matrix[i][k] = 0;
					for(j = k + 1; j < width; ++j) 
					{
						tmp_matrix[i][j] -= quotient * tmp_matrix[k][j];
					}
				}
			}
			min_szie = k;
			for(k = 0; k < min_szie; ++k)
			{
				quotient = tmp_matrix[k][k];
				tmp_matrix[k][k] = 1;
				for(j = k + 1; j < width; ++j) tmp_matrix[k][j] /= quotient;
				for(i = 0; i < k; ++i)
				{
					quotient = tmp_matrix[i][k];
					tmp_matrix[i][k] = 0; 
					for(j = k + 1; j < width; ++j)
					{
						tmp_matrix[i][j] -= quotient * tmp_matrix[k][j];
					}
				}
			}
			k = width - min_szie;
			result.set_size(width, k);
			for(i = 0; i < min_szie; ++i) for(j = 0; j < k; ++j)
			{
				result[index_array[i]][j] = -tmp_matrix[i][min_szie+j];
			}
			for(j = 0; j < k; ++j)
			{
				result[index_array[j + min_szie]][j] = 1;
			}
			return result;
		}

		template <typename value_type>
		value_type det(const matrix <value_type>& var_matrix)
		{
			matrix<value_type> triangle_matrix = trify(var_matrix);
			size_t width = std::min(var_matrix.height(), var_matrix.width());
			value_type result = 1;
			for(size_t i = 0; i < width; ++i)
			{
				result *= triangle_matrix[i][i];
			}
			return result;
		}

		template <typename value_type>
		size_t rank(const matrix <value_type>& var_matrix)
		{
			matrix <value_type> triangle_matrix = trify(var_matrix);
			size_t result = 0;
			size_t width = std::min(var_matrix.height(), var_matrix.width());
			for(; result < width && triangle_matrix[result][result] != 0;  ++result);
			return result;
		}
#pragma endregion

#pragma region _Matrix_Decompose_
		template <typename value_type, typename left_type, typename right_type>
		void trdecom(
			const matrix<value_type>& var_matrix, 
			matrix<left_type>& left_result,
			matrix<right_type>& right_result)
		{	
			size_t i, j, k;
			size_t var_size = var_matrix.height();
			left_result = var_matrix;
			if(var_size != var_matrix.width()) throw exception("The matrix is not square!");
			if(var_size == 0) throw exception("Null matrix!");
			for(k = 0; k < var_size - 1; ++k)
			{
				if(abs(left_result[k][k]) + 1.0 == 1.0)
				{
					throw exception("Triangle decompose failed!");
				}
				for(i = k + 1; i < var_size; ++i)
				{
					left_result[i][k] /= left_result[k][k];
				}
				for(i = k + 1; i < var_size; ++i) for(j = k + 1; j < var_size; ++j)
				{ 
					left_result[i][j] -= left_result[i][k]*left_result[k][j];
				}
			}
			right_result = left_result;
			for(i = 0; i < var_size; ++i)
			{
				for(j = 0; j < i; ++j)
				{
					right_result[i][j] = 0; 
				}
				left_result[i][i] = 1;
				for(j = i + 1; j < var_size; ++j)
				{
					left_result[i][j] = 0; 
				}
			}
		}

		template <typename value_type, typename left_type, typename right_type>
		void qrdecom(
			const matrix<value_type>& var_matrix, 
			matrix<left_type>& left_result,
			matrix<right_type>& right_result)
		{
			size_t i, j, k, l, i_max;				
			size_t height = var_matrix.height(), width = var_matrix.width();
			value_type quotient, conj_quotient, mod;
			ldouble max_abs, tmp_abs;
			array_list<value_type> tmp_vector(std::max(height, width));
			right_result = var_matrix;
			left_result.clear();
			left_result.set_size(height, height);
			for(i = 0; i < height; ++i) left_result[i][i] = 1;
			for(l = 0, k = 0; k < height && l < width; ++l)
			{
				for(max_abs = 0, i = k; i < height; ++i)
				{
					tmp_abs = abs(right_result[i][l]);
					if(tmp_abs > max_abs)
					{
						max_abs = tmp_abs;
						i_max = i;
					}
				}
				if(max_abs + 1.0 == 1.0) continue;
				if(i_max != k)
				{
					for(j = 0; j < width; ++j)
					{
						swap(right_result[k][j], right_result[i_max][j]);
					}
					for(j = 0; j < height; ++j)
					{
						swap(left_result[j][k], left_result[j][i_max]);
					}
				}
				for(i = k + 1; i < height; ++i)
				{
					quotient = right_result[i][l] / right_result[k][l];
					conj_quotient = conj(quotient);
					mod = sqrt(1 + conj_quotient * quotient);
					for(j = l; j < width; ++j)
					{
						tmp_vector[j] = right_result[k][j] + right_result[i][j] * conj_quotient;
						tmp_vector[j] /= mod;
						right_result[i][j] -= quotient * right_result[k][j];
						right_result[i][j] /= mod;
						right_result[k][j] = tmp_vector[j];
					}
					right_result[i][l] = 0; 
					for(j = 0; j < height; ++j)
					{
						tmp_vector[j] = left_result[j][k] + left_result[j][i] * quotient;
						tmp_vector[j] /= mod;
						left_result[j][i] -= conj_quotient * left_result[j][k]
						left_result[j][i] /= mod;
						left_result[j][k] = tmp_vector[j];
					}
				}
				++k;
			}
		}

		template <typename value_type, typename left_type, typename right_type>
		void frdecom(
			const matrix<value_type>& var_matrix, 
			matrix<left_type>& left_result, 
			matrix<right_type>& right_result)
		{
			size_t i, j, k, l, i_max, rank;
			size_t height = var_matrix.height(), width = var_matrix.width();
			ldouble max_abs, tmp_abs;
			value_type quotient;
			right_result = var_matrix;
			left_result.clear();
			left_result.set_size(height, width);
			for(l = 0, k = 0; k < height && l < width; ++l)
			{
				for(max_abs = 0, i = k; i < height; ++i)
				{
					tmp_abs = abs(right_result[i][l]);
					if(tmp_abs > max_abs)
					{
						max_abs = tmp_abs;
						i_max = i;
					}
				}
				if(max_abs + 1.0 == 1.0) continue;
				if(i_max != k)
				{
					for(j = 0; j < width; ++j)
					{
						swap(right_result[k][j], right_result[i_max][j]);
					}
				}
				for(i = k + 1; i < height; ++i)
				{
					quotient = right_result[i][l] / right_result[k][l];
					right_result[i][l] = 0; 
					for(j = l + 1; j < width; ++j)
					{
						right_result[i][j] -= quotient * right_result[k][j];
					}
				}
				++k;
			}
			for(rank = 0, j = 0, l = 0; rank < height; ++rank)
			{
				for(; j < width && right_result[rank][j] == 0;  ++j);
				if(j < width)
				{
					quotient = right_result[rank][j];
					right_result[rank][j] = 1;
					for(k = j + 1; k < width; ++k) 
					{
						right_result[rank][k] /= quotient;
					}
					for(i = 0; i < rank; ++i)
					{
						quotient = right_result[i][j];
						right_result[i][j] = 0; 
						for(k = j + 1; k < width; ++k)
						{ 
							right_result[i][k] -= quotient * right_result[rank][k];
						}
					}
					for(i = 0; i < height; ++i) 
					{
						left_result[i][l] = var_matrix[i][j];						
					}
					 ++l;
				}
				else break;
			}
			left_result.set_size(height, rank);
			right_result.set_size(rank, width);
		}
#pragma endregion
	}
}
#endif
