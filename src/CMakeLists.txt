
SET(MY_LIBRARY_SRCS 
math60082_sor.cpp 
math60082_cn.cpp
math60082_cheb.cpp
math60082_leastSquaresFit.cpp
math60082_gnuplot.cpp
math60082_randomPaths.cpp
math60082_monteCarloLeastSquares.cpp
)

ADD_LIBRARY( math60082 SHARED
${MY_LIBRARY_SRCS}
)

TARGET_LINK_LIBRARIES( math60082 boost_system pthread)

ADD_LIBRARY( math60082_static STATIC
${MY_LIBRARY_SRCS}
)

TARGET_LINK_LIBRARIES( math60082 boost_system pthread)
