#ifndef __GDFMMTRAITS_HPP__
#define __GDFMMTRAITS_HPP__

namespace GDFMM_Traits{

	//Eigen Types
	using MatRow        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using MatCol        = Eigen::MatrixXd;
	using MatUnsRow     = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using MatUnsCol     = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
	using VecUnsCol     = Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>;
	using VecUnsRow     = Eigen::Matrix<unsigned int, 1, Eigen::Dynamic>;
	using VecCol        = Eigen::VectorXd;
	using VecRow        = Eigen::RowVectorXd;

	//Iterators
	using vector_d_iterator 	= std::vector<double>::iterator;
	using vector_d_citerator 	= std::vector<double>::const_iterator;
	using vector_i_iterator 	= std::vector<int>::iterator;
	using vector_i_citerator 	= std::vector<int>::const_iterator;
	using vector_ui_iterator 	= std::vector<unsigned int>::iterator;
	using vector_ui_citerator 	= std::vector<unsigned int>::const_iterator;

}

#endif