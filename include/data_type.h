// Use 2*s for polynomial order to ensure square matrices A and Q.
template <int s>
using VectorP = Eigen::Matrix<double, 2*s, 1>;
template <int s>
using VectorB = Eigen::Matrix<double, 2*s, 1>;
template <int s>
using MatrixA = Eigen::Matrix<double, 2*s, 2*s>;
template <int s>
using MatrixQ = Eigen::Matrix<double, 2*s, 2*s>;
// s-1 is hardcoded as MinimisingPoly is currently only written for fixed position derivatives at each waypoint.
template <int s>
// using MatrixPartition = Eigen::Matrix<double, s-1, s-1>;
using MatrixPartition = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, s, s>;
template <int s>
using BigMatrixPartition = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 2*s, 2*s>;
template <int s>
// using VectorPartition = Eigen::Matrix<double, s-1, 1>;
using VectorPartition = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, s, 1>;
template <int s>
using BigVectorPartition = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 2*s, 1>;

typedef Eigen::Matrix<double, 1, 1> Scalar;

template <int s>
struct MatrixBlock {
	MatrixA<s> A;
	MatrixQ<s> Q;
	MatrixA<s> H;
	VectorB<s> g;
};

template <int s>
struct MatrixSubBCD {
	MatrixPartition<s> B;
	MatrixPartition<s> C;
	MatrixPartition<s> D;
	VectorPartition<s> g0;
	VectorPartition<s> gT;

	VectorB<s> b;
	BigMatrixPartition<s> Z;
	VectorPartition<s> f0;
	VectorPartition<s> fT;
};