#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"

using namespace GDFMM_Traits;

namespace utils{

	//Some constants
	inline constexpr double pi 			     = M_PI;
	inline constexpr double pi_2 			 = M_PI_2;
	inline constexpr double sqrt_2 		     = M_SQRT2;
	inline constexpr double two_over_sqrtpi  = M_2_SQRTPI;
	inline constexpr long double log_2	 	 = M_LN2;
	inline constexpr double sqrt_pi 		 = two_over_sqrtpi*pi_2;
	inline constexpr long double sqrt_2pi 	 = sqrt_pi*sqrt_2;
	inline constexpr long double log_2pi	 = std::log(2*pi);
	inline constexpr long double log_pi	 	 = std::log(pi);

	//------------------------------------------------------------------------------------------------------------------------------------------------------
	//	Sub-Matrix extraction utilities
	//------------------------------------------------------------------------------------------------------------------------------------------------------
	//With the introduction of Eigen 3.4, all these operations can be done much easier and faster. Unfortunately right now RcppEigen supports only Eigen
	// 3.3.7.
	// See https://eigen.tuxfamily.org/dox-devel/group__TutorialSlicingIndexing.html
	
	// The following functor and function are taken from the Eigen manual. https://eigen.tuxfamily.org/dox/TopicCustomizing_NullaryExpr.html.
	// The implement a way to easily select submatrices in Matlab fashion, i.e A([1,2,5],[0,3,6]). 
	// Moreover, removing the const adornal, it is now no longer possible to pass r-value references, i.e SubMatrix({0,1},M) is not possible.
	template<class ArgType, class RowIndexType, class ColIndexType>
	class indexing_functor {
	  const ArgType &m_arg;
	  const RowIndexType &m_rowIndices;
	  const ColIndexType &m_colIndices;
	  public:
	  typedef Eigen::Matrix<typename ArgType::Scalar,
	                 		RowIndexType::SizeAtCompileTime,
	                 		ColIndexType::SizeAtCompileTime,
	                 		ArgType::Flags&Eigen::RowMajorBit?Eigen::RowMajor:Eigen::ColMajor,
	                 		RowIndexType::MaxSizeAtCompileTime,
	                 		ColIndexType::MaxSizeAtCompileTime> MatrixType;
	 
	  indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
	    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
	  {}
	 
	  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
	    return m_arg(m_rowIndices[row], m_colIndices[col]);
	  }
	};
	template <class ArgType, class RowIndexType, class ColIndexType>
	Eigen::CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
	indexing(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
	{
	  typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
	  typedef typename Func::MatrixType MatrixType;
	  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
	}


	//Passing rows to be exctracted and cols to be extracted. May be different and of different EigenType. Can not be Eigen:Array for the moment.
	//Input and Output matrix has to be of the same form, i.e if input is RowMajor, output is RowMajor. Same for ColMajor. Works in both cases.
	//Watch out, it does not check that M is a matrix. It could accept a vector but then there is an error. 
	template<typename Derived_mat, typename Derived_rows, typename Derived_cols>
	Derived_mat SubMatrix(Eigen::MatrixBase<Derived_mat> const & M, Eigen::MatrixBase<Derived_rows> const & idx_rows, Eigen::MatrixBase<Derived_cols> const & idx_cols){

		//Check the type of idx_rows. This makes sure that only eigen vectors are passed or objects that are wrapped into eigen vector by Eigen::Map.
		static_assert( std::is_same_v<Derived_rows, VecUnsCol> || 
					   std::is_same_v<Derived_rows, VecUnsRow> ||
					   std::is_same_v<Derived_rows, Eigen::Map<VecUnsCol>> ||
					   std::is_same_v<Derived_rows, Eigen::Map<const VecUnsCol>> ||
					   std::is_same_v<Derived_rows, Eigen::Map<VecUnsRow>> ||
					   std::is_same_v<Derived_rows, Eigen::Map<const VecUnsRow>> ,
					  "______ ERROR, invalid Input Type requested in SubMatrix. Can handle only VecUnsRow, VecUnsCol  _____");

		//Check the type of idx_cols. This makes sure that only eigen vectors are passed or objects that are wrapped into eigen vector by Eigen::Map.
		static_assert( std::is_same_v<Derived_cols, VecUnsCol> || 
					   std::is_same_v<Derived_cols, VecUnsRow> ||
					   std::is_same_v<Derived_cols, Eigen::Map<VecUnsCol>> ||
					   std::is_same_v<Derived_cols, Eigen::Map<const VecUnsCol>> ||
					   std::is_same_v<Derived_cols, Eigen::Map<VecUnsRow>> ||
					   std::is_same_v<Derived_cols, Eigen::Map<const VecUnsRow>> ,
					  "______ ERROR, invalid Input Type requested in SubMatrix. Can handle only VecUnsRow, VecCol  _____");

		if(idx_rows.size() > M.rows() || idx_cols.size() > M.cols()) //Check dimension
			throw std::runtime_error("Indeces exceed matrix dimension");

		Derived_mat res = indexing(M, idx_rows, idx_cols);
		return res;
	} 
	
	//Function overload for different specifications
	template<typename Derived_mat, typename Derived_cols>
	Derived_mat SubMatrix(Eigen::MatrixBase<Derived_mat> const & M, unsigned int const & idx_rows, Eigen::MatrixBase<Derived_cols> const & idx_cols){
		VecUnsRow eigen_idx_row(VecUnsRow::Constant(1, idx_rows));
		return SubMatrix(M, eigen_idx_row, idx_cols);
	}

	template<typename Derived_mat, typename Derived_rows>
	Derived_mat SubMatrix(Eigen::MatrixBase<Derived_mat> const & M, Eigen::MatrixBase<Derived_rows> const & idx_rows, unsigned int const & idx_cols){
		return SubMatrix(M, idx_rows, Eigen::Map<const VecUnsCol> (&idx_cols, 1));
	}
}

#endif