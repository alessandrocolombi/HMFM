#ifndef __GDFMMTRAITS_HPP__
#define __GDFMMTRAITS_HPP__

namespace GDFMM_Traits{

	//Eigen Types
	using MatRow        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using MatCol        = Eigen::MatrixXd;
	using MatUnsRow     = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using MatUnsCol     = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
	using VecUnsCol     = Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>;
	using VecCol        = Eigen::VectorXd;
	using VecRow        = Eigen::RowVectorXd;

}

#endif