#ifndef __TESTEIGEN_HPP__
#define __TESTEIGEN_HPP__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"
#include "utils.h"

using namespace GDFMM_Traits;

//' Eigen library usage example
//'
//' @export
// [[Rcpp::export]]
void EigenTest(){
	int p = 4;

	// Eigen type definition and initialization
	Rcpp::Rcout<<"Some matrices "<<std::endl;
	MatCol M1; //Define but does not initialize. Avoid it
	MatCol M2 = MatCol::Constant(p,p,1.0); 	//Define and initialize a pxp matrix of all ones.
	MatCol M3( MatCol::Constant(p,p,1.0) ); //Same as before but more efficient
	MatRow I1(MatRow::Identity(p,p)); 		//RowMajor pxp identity matrix
	MatCol I2(MatCol::Identity(p,p)); 		//ColMajor pxp identity matrix
	MatCol M4(MatCol::Random(p,p)); 		//ColMajor pxp matrix filled with random numbers
	M1.setRandom(p+p,p+p); 					//resize M1 and fill it with random numbers

	//Stream operator is available for eigen objs
	Rcpp::Rcout<<"Print Identity: "<<std::endl<<I1<<std::endl;
	Rcpp::Rcout<<"Set the value of the diagonal elements"<<std::endl;
	Rcpp::Rcout<<"M3.diagonal().array() = 5.5"<<std::endl;
	M3.diagonal().array() = 5.5;

	Rcpp::Rcout<<"Some vectors"<<std::endl;
	VecCol ones(VecCol::Constant(p,1.0));
	VecRow zeros(VecRow::Constant(p,1.0));
	VecRow zeros2(VecRow::Zero(p)); //as before
	VecRow v1(VecRow::Random(p));
	VecCol v2(VecCol::Random(p));
	VecUnsRow v3(VecUnsRow::Constant(p,5)); //vectors of unsigned intergers
	VecUnsCol v4(VecUnsCol::Constant(p,5));	//vectors of unsigned intergers

	
	Rcpp::Rcout<<"Resize and fill by hand"<<std::endl;	
	Rcpp::Rcout<<"v3 = "<<v3<<std::endl;
	Rcpp::Rcout<<"v4 = "<<v4<<std::endl;
	v3.resize(2); //resize the vector and delete olds elements
	v4.conservativeResize(2); //resize the vector and keeps olds elements (less efficient)
	Rcpp::Rcout<<"v3 = "<<v3<<std::endl; //values are trash
	Rcpp::Rcout<<"v4 = "<<v4<<std::endl;
	//Fill vectors by hands (works also for matrices)
	v3<<0,1; 
	v4<<0,2;
	Rcpp::Rcout<<"v3 = "<<v3<<std::endl;
	Rcpp::Rcout<<"v4 = "<<v4<<std::endl;

	Rcpp::Rcout<<"Get the diagonal of a matrix and treat it as a vector: A.diagonal() --> works like diag(A) in R"<<std::endl;
	Rcpp::Rcout<<"ones + I1.diagonal() = :"<<std::endl<<ones + I1.diagonal() <<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;

	Rcpp::Rcout<<"##########################"<<std::endl;
	Rcpp::Rcout<<"    Cast Eigen matrix    " <<std::endl;
	Rcpp::Rcout<<"##########################"<<std::endl;
	Rcpp::Rcout<<"It is possible to cast eigen objects into a new type."<<std::endl;
	Rcpp::Rcout<<"Cast a matrix of int into double"<<std::endl;
	MatIntRow Aint(MatIntRow::Constant(p,p,5));
	Rcpp::Rcout<<"I1*Aint.cast<double>():"<<std::endl<<I1*Aint.cast<double>()<<std::endl;

	Rcpp::Rcout<<std::endl<<std::endl;
	Rcpp::Rcout<<"##########################"<<std::endl;
	Rcpp::Rcout<<"     Basic operations     " <<std::endl;
	Rcpp::Rcout<<"##########################"<<std::endl;

	Rcpp::Rcout<<"---- Access elements ----"<<std::endl;
	Rcpp::Rcout<<"Use call operator() to access elements"<<std::endl;
	std::cout<<"M4(1,1):"<<std::endl<<M4(1,1)<<std::endl;
	Rcpp::Rcout<<"ones(1):"<<std::endl<<ones(1)<<std::endl;
	//get dimensions for matrix
	Rcpp::Rcout<<"---- get dimensions ----"<<std::endl;
	Rcpp::Rcout<<"Example for matrices"<<std::endl;
	int n_elem = M4.size(); //the number of elements, i.e p*p
	int n_rows = M4.rows(); //the number of rows, i.e p
	int n_cols = M4.cols(); //the number of cols, i.e p
	Rcpp::Rcout<<"n_elem"<<n_elem<<std::endl;
	Rcpp::Rcout<<"n_rows"<<n_rows<<std::endl;
	Rcpp::Rcout<<"n_cols"<<n_cols<<std::endl;
	//get dimension for vectors, internally they are matrices, so it works the same
	Rcpp::Rcout<<"Example for vectors"<<std::endl;
	n_elem = v1.size(); //the number of elements, i.e p*p
	n_rows = v1.rows(); //the number of rows, i.e p
	n_cols = v1.cols(); //the number of cols, i.e p
	Rcpp::Rcout<<"n_elem"<<n_elem<<std::endl;
	Rcpp::Rcout<<"n_rows"<<n_rows<<std::endl;
	Rcpp::Rcout<<"n_cols"<<n_cols<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;

	Rcpp::Rcout<<"---- Transpose operation ----"<<std::endl;
	MatCol M5 = M4.transpose(); //ok
	Rcpp::Rcout<<"M4:"<<std::endl<<M4<<std::endl;
	Rcpp::Rcout<<"M5 = M4.transpose():"<<std::endl<<M5<<std::endl;
	Rcpp::Rcout<<"Some examples"<<std::endl;
	Rcpp::Rcout<<"M5 = M5 + M4.transpose() is ok"<<std::endl;
	M5 = M5 + M4.transpose(); // ok
	Rcpp::Rcout<<"M5:"<<std::endl<<M5<<std::endl;
	Rcpp::Rcout<<"M5 = M4 + M5.transpose() can not be done"<<std::endl;
	//M5 = M4 + M5.transpose(); //NO!!!
	Rcpp::Rcout<<"M5:"<<std::endl<<M5<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;

	Rcpp::Rcout<<"---- Linear Algebra operations ----"<<std::endl;
	Rcpp::Rcout<<"Can mix RowMajor and ColMajor. It is possible but may be inefficient, choose carefully the storage method"<<std::endl;
	MatCol res = 3*M1 + 4*I1; 
	Rcpp::Rcout<<"res = 3*M1 + 4*I1 :"<<std::endl<<res<<std::endl;
	Rcpp::Rcout<<"Products"<<std::endl;
	MatCol res2 = M2*M4;
	Rcpp::Rcout<<"res2 = M2*M4 :"<<std::endl<<res2<<std::endl;
	VecCol res_vec = M2*ones;
	Rcpp::Rcout<<"res_vec = M2*ones :"<<std::endl<<res_vec<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;

	Rcpp::Rcout<<"----Eigen Reductions, from Eigen object to scalar ----"<<std::endl<<"See https://eigen.tuxfamily.org/dox/group__TutorialReductionsVisitorsBroadcasting.html for more"<<std::endl;
	Rcpp::Rcout<<"They works for matrices and vectors the same way"<<std::endl;
	Rcpp::Rcout<<"sum"<<std::endl;
	Rcpp::Rcout<<"M5.sum() = :"<<std::endl<<M5.sum()<<std::endl;
	Rcpp::Rcout<<"v1.sum() =:"<<std::endl<<v1.sum()<<std::endl;
	Rcpp::Rcout<<"product of elements"<<std::endl;
	Rcpp::Rcout<<"v1.prod() =:"<<std::endl<<v1.prod()<<std::endl;
	Rcpp::Rcout<<"mean "<<std::endl;
	Rcpp::Rcout<<"v1.mean() =:"<<std::endl<<v1.mean()<<std::endl;
	Rcpp::Rcout<<"minimum"<<std::endl;
	Rcpp::Rcout<<"v1.minCoeff() =:"<<std::endl<<v1.minCoeff()<<std::endl;
	Rcpp::Rcout<<"maximum"<<std::endl;
	Rcpp::Rcout<<"v1.maxCoeff() =:"<<std::endl<<v1.maxCoeff()<<std::endl;
	Rcpp::Rcout<<"Trace (for matrices)"<<std::endl;
	Rcpp::Rcout<<"M5.trace() = :"<<std::endl<<M5.trace()<<std::endl;
	Rcpp::Rcout<<"Inner product between vectors"<<std::endl;
	Rcpp::Rcout<<"ones:"<<std::endl<<ones<<std::endl;
	Rcpp::Rcout<<"v2:"<<std::endl<<v2<<std::endl;
	Rcpp::Rcout<<"v2.dot(ones):"<<std::endl<<v2.dot(ones)<<std::endl;
	Rcpp::Rcout<<"ones.dot(v2):"<<std::endl<<ones.dot(v2)<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;


	Rcpp::Rcout<<"---- Norms and absolute value ----"<<std::endl;
	Rcpp::Rcout<<"M5:"<<std::endl<<M5<<std::endl;
	Rcpp::Rcout<<"Norm 1 = max_j sum_i (|a_ij|)"<<std::endl;
	Rcpp::Rcout<<"M5.lpNorm<1>():"<<std::endl<<M5.lpNorm<1>()<<std::endl;
	Rcpp::Rcout<<"Norm infinity = max_i sum_j (|a_ij|)"<<std::endl;
	Rcpp::Rcout<<"M5.lpNorm<Eigen::Infinity>():"<<std::endl<<M5.lpNorm<Eigen::Infinity>()<<std::endl;
	Rcpp::Rcout<<"Norm 2 = sum_ij (|a_ij|^2)"<<std::endl;
	Rcpp::Rcout<<"M5.squaredNorm():"<<std::endl<<M5.squaredNorm()<<std::endl;
	VecCol v_negative = VecCol::Constant(p,-2.2);
	Rcpp::Rcout<<"v_negative:"<<std::endl<<v_negative<<std::endl;
	VecCol v_positive = v_negative.cwiseAbs();
	Rcpp::Rcout<<"v_positive = v_negative.cwiseAbs():"<<std::endl<<v_positive<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;

	Rcpp::Rcout<<"---- Eigen PartialReduction, from Matrix to vector. ----"<<std::endl<<"They are operations done row by row or columns by columns"<<std::endl;
	Rcpp::Rcout<<".colwise() means that the operation has to be done in each column. Same for .rowwise()"<<std::endl;
	Rcpp::Rcout<<"M5:"<<std::endl<<M5<<std::endl;
	Rcpp::Rcout<<"M5.rowwise().maxCoeff():"<<std::endl<<M5.rowwise().maxCoeff()<<std::endl;
	Rcpp::Rcout<<"M5.colwise().maxCoeff():"<<std::endl<<M5.colwise().maxCoeff()<<std::endl;
	Rcpp::Rcout<<"M5.rowwise().sum():"<<std::endl<<M5.rowwise().sum()<<std::endl;
	Rcpp::Rcout<<"M5.colwise().sum():"<<std::endl<<M5.colwise().sum()<<std::endl;
	Rcpp::Rcout<<"Broadcast: sum a vector to each column or row. Watch out, if I use .colwise() I need to sum a column vector. I use .rowwise() I need to sum a row vector."<<std::endl;
	Rcpp::Rcout<<"M5:"<<std::endl<<M5<<std::endl;
	Rcpp::Rcout<<"ones:"<<std::endl<<ones<<std::endl;
	M5.colwise()+=ones;
	Rcpp::Rcout<<"M5.colwise()+=ones.transpose():"<<std::endl<<M5<<std::endl;
	M5.rowwise()+=ones.transpose();
	Rcpp::Rcout<<"M5.rowwise()+=ones.transpose():"<<std::endl<<M5<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;

	Rcpp::Rcout<<"##########################"<<std::endl;
	Rcpp::Rcout<<"  Sub-Matrix extraction   " <<std::endl;
	Rcpp::Rcout<<"##########################"<<std::endl;

	Rcpp::Rcout<<"Get column j, get row j: use .col(j) and .row(j) methods"<<std::endl;
	Rcpp::Rcout<<"I1.col(2):"<<std::endl<<I1.col(2)<<std::endl;
	Rcpp::Rcout<<"I1.row(2):"<<std::endl<<I1.row(2)<<std::endl;
	Rcpp::Rcout<<"General function for matrices, like A[rows,cols] in R. Does not work for vector, v[idx] is not implemented for idx generic."<<std::endl;
	VecUnsRow idx_rows(v3); //vectors of unsigned intergers
	VecUnsCol idx_cols(v4);	//vectors of unsigned intergers

	Rcpp::Rcout<<"idx_rows = "<<idx_rows<<std::endl;
	Rcpp::Rcout<<"idx_cols = "<<idx_cols<<std::endl;
	Rcpp::Rcout<<"SubMatrix extraction (RowMajor input and output):"<<std::endl<<"I = "<<std::endl<<I1<<std::endl;
	//This is like writing I1[idx_rows,idx_cols] in R.
	MatRow SubI = utils::SubMatrix(I1,idx_rows,idx_cols);
	Rcpp::Rcout<<"SubI = "<<std::endl<<SubI<<std::endl;
	//MatRow SubI_2 = utils::SubMatrix(I1,{0,1},idx_cols); //this does not work
	MatRow SubI_2 = utils::SubMatrix(I1,2,idx_cols); 
	Rcpp::Rcout<<"SubI_2 = I1[2,idx_cols] = "<<std::endl<<SubI_2<<std::endl;
	MatCol SubI_3 = utils::SubMatrix(I1,idx_rows,1); 
	Rcpp::Rcout<<"SubI_3 = I1[idx_rows,1] = "<<std::endl<<SubI_3<<std::endl;

	Rcpp::Rcout<<"SubMatrix extraction (ColMajor input and output):"<<std::endl<<"M4 = "<<std::endl<<M4<<std::endl;
	MatCol SubM4 = utils::SubMatrix(M4,idx_rows,idx_cols);
	Rcpp::Rcout<<"SubM4 = "<<std::endl<<SubM4<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;

	Rcpp::Rcout<<"Works also for vectors, just set equal to 0 one of the indicies."<<std::endl;
	VecCol r1(VecCol::Random(p));
	std::cout<<"r1:"<<std::endl<<r1<<std::endl;
	std::cout<<"utils::SubMatrix(r1,idx_rows,0):"<<std::endl<<utils::SubMatrix(r1,idx_rows,0)<<std::endl;
	VecRow r2(VecRow::Random(p));
	Rcpp::Rcout<<"r2:"<<std::endl<<r2<<std::endl;
	std::cout<<"utils::SubMatrix(r2,0,idx_cols):"<<std::endl<<utils::SubMatrix(r2,0,idx_cols)<<std::endl;
	

	Rcpp::Rcout<<"---- Block Operations for matrices (SubMatrix is better) ----"<<std::endl;
	Rcpp::Rcout<<"A.block(i,j,dim_row,dim_col) -> gets a pxq sub_matrix starting from element (i,j), which is the top left element. Examples:"<<std::endl;
	Rcpp::Rcout<<"M4.block(1,0,3,2):"<<std::endl<<M4.block(1,0,3,2)<<std::endl;
	Rcpp::Rcout<<"There are shortcuts if the block starts in corner of the matrix. Top/Bottom and Left/Right corner"<<std::endl;
	Rcpp::Rcout<<"M4.topRightCorner(3,2):"<<std::endl<<M4.topRightCorner(3,2)<<std::endl;
	Rcpp::Rcout<<"M4.topLeftCorner(3,2):"<<std::endl<<M4.topLeftCorner(3,2)<<std::endl;
	Rcpp::Rcout<<"M4.bottomRightCorner(3,2):"<<std::endl<<M4.bottomRightCorner(3,2)<<std::endl;
	Rcpp::Rcout<<"M4.bottomLeftCorner(3,2):"<<std::endl<<M4.bottomLeftCorner(3,2)<<std::endl;
	Rcpp::Rcout<<"Additional shortcuts to get the first/last rows/columns"<<std::endl;
	Rcpp::Rcout<<"M5.topRows(2):"<<std::endl<<M5.topRows(2)<<std::endl;
	Rcpp::Rcout<<"M5.bottomRows(2):"<<std::endl<<M5.bottomRows(2)<<std::endl;
	Rcpp::Rcout<<"M5.leftCols(2):"<<std::endl<<M5.leftCols(2)<<std::endl;
	Rcpp::Rcout<<"M5.rightCols(2):"<<std::endl<<M5.rightCols(2)<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;

	Rcpp::Rcout<<"---- Block Operations for vectors ----"<<std::endl;
	VecRow v6(VecRow::Random(p));
	std::cout<<"v6:"<<std::endl<<v6<<std::endl;
	Rcpp::Rcout<<"v.segment(i,n) -> gets n elements starting from element i"<<std::endl;
	Rcpp::Rcout<<"v6.segment(1,2):"<<std::endl<<v6.segment(1,2)<<std::endl;
	Rcpp::Rcout<<"Shortcuts to the the first/last n elements"<<std::endl;
	Rcpp::Rcout<<"v6.head(2):"<<std::endl<<v6.head(2)<<std::endl;
	Rcpp::Rcout<<"v6.tail(2):"<<std::endl<<v6.tail(2)<<std::endl;
	Rcpp::Rcout<<std::endl<<std::endl;
	
	Rcpp::Rcout<<"##########################"<<std::endl;
	Rcpp::Rcout<<"     Matrix Inversion     " <<std::endl;
	Rcpp::Rcout<<"##########################"<<std::endl;

	MatRow Irow(MatRow::Identity(p,p));
	MatRow A(MatRow::Random(p,p));
	Rcpp::Rcout<<"A:"<<std::endl<<A<<std::endl;
	MatRow InvA = A.inverse(); //uses LU factorization 
	Rcpp::Rcout<<"InvA:"<<std::endl<<InvA<<std::endl;
	Rcpp::Rcout<<"If A is symmetric positive definite, it is better to use the Cholesky decomposition. More efficient than LU"<<std::endl;
	/*
		InvA = A.llt().solve(Irow); //Computed the Cholesky decomposition of A and solves a linear system. Does not change much if I ir RowMajor or ColMajor and if A is Row/ColMajor
		Rcpp::Rcout<<"InvA:"<<std::endl<<InvA<<std::endl;
	*/

	Rcpp::Rcout<<std::endl<<std::endl;
	
	Rcpp::Rcout<<"##########################"<<std::endl;
	Rcpp::Rcout<<"     Useful Example     " <<std::endl;
	Rcpp::Rcout<<"##########################"<<std::endl;
	int r(5);
	int n(10);
	MatCol data(MatCol::Random(r,n)); //rxn matrix. Think of data as a matrix of n r-dimensional variables. y_1,...,y_n where y_i is in R^r
	Rcpp::Rcout<<"data:"<<std::endl<<data<<std::endl;
	Rcpp::Rcout<<"C = A.cwiseProduct(B) does the component wise product of matricies A and B, i.e c_ij = a_ij*b_ij"<<std::endl;
	//std::cout<<"data.cwiseProduct(data):"<<std::endl<<data.cwiseProduct(data)<<std::endl;
	Rcpp::Rcout<<"It is useful to combine operations to get efficient results without saving any partial results. "<<std::endl;
	Rcpp::Rcout<<"For example: data.cwiseProduct(data).sum() = "<<data.cwiseProduct(data).sum()<<std::endl<< "is the sum of all the elements squared that is the same of computing sum_i( t(data[,i])*data[,i] ) = sum_i(<y_i,y_i>), where <,> is the inner product ."<<std::endl;
	Rcpp::Rcout<<"Indeed, compute the same with the for loop.."<<std::endl;
	double S = 0.0;
	for(int i = 0; i < n; ++i){
		S += data.col(i).dot(data.col(i));
	}
	Rcpp::Rcout<<"S = sum_i(data.col(i).dot(data.col(i))):"<<std::endl<<S<<std::endl;
	Rcpp::Rcout<<"Avoid for loop and use eigen operations!"<<std::endl;
}


#endif