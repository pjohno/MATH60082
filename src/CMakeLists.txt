
SET(MY_LIBRARY_SRCS 
math60082_sor.cpp 
math60082_cn.cpp
math60082_cheb.cpp
math60082_leastSquaresFit.cpp
math60082_gnuplot.cpp
math60082_randomPaths.cpp
math60082_gsl_fit.cpp
math60082_monteCarloLeastSquares.cpp
math60082_tridag.cpp
math60082_lagrange_interp.cpp
math60082_stencils.cpp
math60082_blackScholes.cpp
math60082_parabollicPDE.cpp
math60082_parabollicPDE_examples.cpp
math60082_fokker_planck.cpp
math60082_mVector.cpp
)

ADD_LIBRARY( math60082 SHARED
${MY_LIBRARY_SRCS}
)

TARGET_LINK_LIBRARIES( math60082 boost_system pthread gsl cblas)

ADD_LIBRARY( math60082_static STATIC
${MY_LIBRARY_SRCS}
)

TARGET_LINK_LIBRARIES( math60082_static boost_system pthread gsl cblas)
