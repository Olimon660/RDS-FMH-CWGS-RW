#include <cmath>
#include <vector>
#include <Eigen/Dense>

namespace DreamLab {
	
    struct NonInveribleMatix : std::runtime_error
    {
        NonInveribleMatix() : std::runtime_error("NonInveribleMatix") {}
    };
    
    using Matrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;
    
    typedef int Seed;
    
    // Beta coefficient
    typedef double Beta;
    
    // Number of iterations
    typedef unsigned Convergence;
    
    // Solution converged?
    typedef bool Converged;
    
    typedef std::vector<int> Permutation;
    
    // Regression Statistics
    typedef std::tuple<Matrix, Convergence, Converged> RStats;
    
    // Genotype statistics
    typedef std::tuple<Beta, Convergence, Converged> GStats;
    
    class Generator {
    public:
        Generator(Seed s1 = 3141, Seed s2 = 5926, Seed s3 = 5358) : _s1(s1), _s2(s2) , _s3(s3) {}        
        Generator(const Generator&) = default;
        
        inline double runif()
        {
            _s1 = (171 * _s1) % 30269;
            _s2 = (172 * _s2) % 30307;
            _s3 = (170 * _s3) % 30323;
            double tmp;
            return modf(_s1 / 30269. + _s2 / 30307. + _s3 / 30323., &tmp);
        }
        
        inline int rint(int lower, int upper)
        {
            return lower + int(runif() * (upper - lower));
        }
        
        inline Permutation rperm(unsigned size)
        {
            Permutation x;
            x.resize(size);
            
            for (auto i = 0; i < x.size(); i++) {
                x[i] = i;
            }
            
            for (auto i = 0; i < size - 1; i++) {
                const auto j = rint(i, size);
                
                if (i != j) {
                    std::swap(x[i], x[j]);
                }
            }
            
            return x;
        }
        
        inline Seed s1() const { return _s1; }
        inline Seed s2() const { return _s2; }
        inline Seed s3() const { return _s3; }
        
    private:
        int _s1, _s2, _s3;
    };
    
    struct Linear
    {
        template <class T>
        static Eigen::Matrix<typename T::Scalar, T::ColsAtCompileTime, T::RowsAtCompileTime>
        pseudoinverse(const T& mat, typename T::Scalar tolerance = typename T::Scalar{ 1e-4 })
        {
            typedef typename T::Scalar Scalar;
            auto svd = mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
            const auto& singularValues = svd.singularValues();
            Eigen::Matrix<Scalar, T::ColsAtCompileTime, T::RowsAtCompileTime> singularValuesInv(mat.cols(), mat.rows());
            singularValuesInv.setZero();
            for (auto i = 0; i < singularValues.size(); i++) {
                if (singularValues(i) > tolerance) {
                    singularValuesInv(i, i) = Scalar{ 1 } / singularValues(i);
                }
                else {
                    singularValuesInv(i, i) = Scalar{ 0 };
                }
            }
            return svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
        }
        
        template <class T>
        static RStats linear(const T& Y, const T& X)
        {
            if (Y.cols() > 1) {
                throw std::runtime_error("Invalid dependent variable. It must have a column size of 1.");
            }
            else if (X.rows() != Y.rows()) {
                throw std::runtime_error("Invalid inputs. Variables must have identical number of rows.");
            }
            
            return RStats(Linear::pseudoinverse(X) * Y, 1, true);
        }
        
        template <class T>
        static GStats solve(const T& Y, const T& X)
        {
            const auto r = Linear::linear(Y, X);
            return GStats(std::get<0>(r)(1), std::get<1>(r), std::get<2>(r));
        }
    };
    
    struct Logistic
    {
        static RStats logistic(const Matrix& Y, const Matrix& X, long double min_delta = 1e-6, int max_iters = 100)
        {
            if (Y.cols() > 1) {
                throw std::runtime_error("Invalid dependent variable. It must have a column size of 1.");
            }
            else if (X.rows() != Y.rows()) {
                throw std::runtime_error("Invalid inputs. Variables must have identical number of rows.");
            }
            else if ((Y.array() != 0 && Y.array() != 1).any()) {
                throw std::runtime_error("Invalid dependent variable. The values must be 1 or 0.");
            }
            
            Matrix beta;
            beta.setZero(X.cols(), 1);
            Matrix ones;
            ones.setOnes(X.rows(), 1);
            
            Matrix mu = X * beta;
            
            auto logis = [](long double x) -> long double { if (x < -700) x = -700;
                if (x > 700) x = 700;
                return 1. / (1. + expl(-x)); };
            Matrix pi = mu.unaryExpr(logis);
            
            auto logexp = [](long double x) -> long double { return logl(1. + expl(-x)); };
            long double loglik = -(mu.unaryExpr(logexp) + mu.cwiseProduct(ones - Y)).sum();
            
            long double stepsize = 1.;
            
            int i = 0;
            for (i = 0; i < max_iters; ++i) {
                Matrix Xa = X.cwiseProduct(pi.unaryExpr([](long double x) { return x * (1 - x); }) * Matrix().setOnes(1, X.cols()));
                Matrix delta = (X.transpose() * Xa).inverse() * X.transpose() * (Y - pi);
                Matrix new_beta;
                Matrix new_mu;
                Matrix new_pi;
                long double new_loglik;
                
                while (true) {
                    new_beta = beta + delta * stepsize;
                    new_mu = X * new_beta;
                    new_pi = new_mu.unaryExpr(logis);
                    new_loglik = -((new_mu.unaryExpr(logexp) + new_mu.cwiseProduct(ones - Y)).sum());
                    
                    if (new_loglik > loglik) {
                        if (new_loglik - loglik < min_delta) {
                            return RStats(new_beta, i + 1, true);
                        }
                        
                        break;
                    }
                    else {
                        if (loglik - new_loglik < min_delta) {
                            return RStats(beta, i + 1, true);
                        }
                    }
                    stepsize /= 2.;
                }
                pi = new_pi;
                beta = new_beta;
                loglik = new_loglik;
            }
            
            return RStats(beta, i + 1, false);
        }
        
        template <class T>
        static GStats solve(const T& Y, const T& X)
        {
            const auto r = Logistic::logistic(Y, X);
            return GStats(std::get<0>(r)(1, 0), std::get<1>(r), std::get<2>(r));
        }
    };
    
    struct Firth
    {
        static long double LD(const Matrix& m)
        {
            Eigen::PartialPivLU<Matrix> lu(m);
            auto& LU = lu.matrixLU();
            long double ld = 0.;
            long double c = lu.permutationP().determinant(); // sign
            for (int i = 0; i < LU.rows(); ++i) {
                long double ith_diag = LU(i, i);
                if (ith_diag < 0.)
                    c *= -1;
                ld += logl(abs(ith_diag));
            }
            ld += logl(c);
            return ld;
        }
        
        static RStats firth(const Matrix& Y, const Matrix& X, long double min_delta = 1e-8, int max_iters = 100)
        {
            if (Y.cols() > 1) {
                throw std::runtime_error("Invalid dependent variable. It must have a column size of 1.");
            }
            else if (X.rows() != Y.rows()) {
                throw std::runtime_error("Invalid inputs. Variables must have identical number of rows.");
            }
            else if ((Y.array() != 0 && Y.array() != 1).any()) {
                throw std::runtime_error("Invalid dependent variable. The values must be 1 or 0.");
            }
            
            int maxstep = 5;
            int maxhs = 5;
            
            Matrix beta;
            beta.setZero(X.cols(), 1);
            Matrix ones;
            ones.setOnes(X.rows(), 1);
            
            auto logis = [](long double x) -> long double { if (x < -700) x = -700;
                if (x > 700) x = 700;
                return 1. / (1. + expl(-x)); };
            auto logexp = [](long double x) -> long double { return logl(1. + expl(-x)); };
            
            Matrix Xbeta = X * beta;
            Matrix pi = Xbeta.unaryExpr(logis);
            long double loglik = -(Xbeta.unaryExpr(logexp) + Xbeta.cwiseProduct(ones - Y)).sum();
            
            Matrix Xw2 = X.cwiseProduct(pi.unaryExpr([](long double x) { return sqrtl(x * (1 - x)); }) * Matrix().setOnes(1, X.cols()));
            Matrix Fisher = Xw2.transpose() * Xw2;
            Matrix covs = Fisher.inverse();
            
            long double logdet = LD(Fisher); /*logl(Fisher.determinant());*/
            
            /*
             * Non-invertible matrix if one or more columns are perfectly linearly dependent
             */
            
            if (isnan(logdet) || isinf(logdet))
            {
				throw NonInveribleMatix();
            }
            
            loglik += logdet / 2.;
            int i = 0;
            for (i = 0; i < max_iters; ++i) {
                long double loglik_old = loglik;
                Matrix beta_old = beta;
                
                Xw2 = X.cwiseProduct(pi.unaryExpr([](long double x) { return sqrtl(x * (1 - x)); }) * Matrix().setOnes(1, X.cols()));
                
                Fisher = Xw2.transpose() * Xw2;
                covs = Fisher.inverse();
                Matrix Hdiag = (Xw2 * covs * Xw2.transpose()).diagonal();
                
                Matrix residp = Y - pi + Hdiag.cwiseProduct((0.5 * ones) - pi);
                Matrix U_star = (X.transpose() * residp).transpose();
                Matrix XX_XW2 = (X.transpose() * pi.unaryExpr([](long double x) { return sqrtl(x * (1. - x)); }).col(0).asDiagonal());
                
                Matrix XX_Fischer = XX_XW2 * XX_XW2.transpose();
                Matrix XX_covs = XX_Fischer.inverse();
                Matrix delta = XX_covs * U_star.transpose();
                long double loglik_change = 0;
                
                long double mx = delta.cwiseAbs().maxCoeff() / maxstep;
                if (mx > 1.)
                    delta = delta / mx;
                
                double stepsize = 1.;
                beta += delta;
                
                for (int j = 0; j < maxhs; ++j) {
                    Xbeta = X * beta;
                    pi = Xbeta.unaryExpr(logis);
                    
                    loglik = -((Xbeta.unaryExpr(logexp) + Xbeta.cwiseProduct(ones - Y)).sum());
                    Xw2 = X.cwiseProduct(pi.unaryExpr([](long double x) { return sqrtl(x * (1 - x)); }) * Matrix().setOnes(1, X.cols()));
                    
                    Fisher = Xw2.transpose() * Xw2;
                    logdet = LD(Fisher); /*logl(Fisher.determinant());*/
                    
                    loglik += logdet / 2.;
                    
                    loglik_change = loglik - loglik_old;
                    
                    if (loglik >= loglik_old)
                        break;
                    
                    stepsize /= 2.;
                    beta -= delta / stepsize;
                }
                if (loglik_change < min_delta)
                    return RStats(beta, i + 1, true);
            }
            return RStats(beta, i + 1, false);
        }
        
        template <class T> static GStats solve(const T& Y, const T& X)
        {
            const auto r = Firth::firth(Y, X);
            return GStats(std::get<0>(r)(1, 0), std::get<1>(r), std::get<2>(r));
        }
    };
    
    /*
     * Fit model (initialised with values of y and X), and return the fitted genotype coefficient, beta[1].
     *
     *   model: An initialised model of class Linear, Logistic, or Firth.
     *   kwargs: keyword arguments to the model's solve method.
     *
     * Returns (beta[1], conv) where beta[1] is a floating point number, and the boolean conv is True if the
     * fit converged, else false.
     */
    template <typename Model> GStats _gt_statistic(const Matrix& Y, const Matrix& X)
    {
        try
        {
            return Model().solve(Y, X);
        }
        catch (const NonInveribleMatix &)
        {
            return GStats(NAN, 0, false);
        }
    }
    
    /*
     * Calculate the genotype coefficient for model on permuted data.
     *
     *   model: An initialised model of class Linear, Logistic, or Firth.
     *   prng: A seeded prng.
     *   kwargs: keyword arguments to the model's solve method.
     *
     * Returns beta[1] if the fit converged, else None.
     */
    
    template <typename Model> GStats _gt_permutation(const Matrix& Y, const Matrix& X, Generator& g)
    {
        auto Y_ = Y;
        
        auto perm = [&](int size, Generator& g) {
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(size);
            perm.setIdentity();
            for (int i = 0; i < size; ++i) {
                auto j = g.rint(i, size);
                std::swap(perm.indices()[i], perm.indices()[j]);
            }
            return perm;
        };
        
        // Apply permuation
        Y_ = perm(Y.rows(), g).transpose() * Y;
        
        try
        {
            return Model().solve(Y_, X);
        }
        catch (const NonInveribleMatix &)
        {
            return GStats(NAN, 0, false);
        }
    }
    
    /*
     * Sample the null distribution of the genotype coefficient for model by permutation of sample labels.
     *
     *   model: An initialised model of class Linear, Logistic, or Firth.
     *   B: number of permutations
     *   kwargs: keyword arguments to the model's solve method.
     *
     * Returns a list (with B elements) of the permuted coefficients. Elements of this list may be None,
     * if the corresponding model fit failed.
     */
    
    template <typename Model>
    std::vector<GStats> _gt_statistic_null(const Matrix& Y, const Matrix& X, std::size_t B, int s1, int s2, int s3)
    {
        Generator g1(s1, s2, s3);
        
        std::vector<std::vector<int> > seeds(3, std::vector<int>(B));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < B; ++j)
                seeds[i][j] = g1.rint(1, 30001);
        
        std::vector<GStats> x;
        for (auto i = 0; i < B; i++) {
            Generator g2(seeds[0][i], seeds[1][i], seeds[2][i]);
            x.push_back(_gt_permutation<Model>(Y, X, g2));
        }
        return x;
    }
}
