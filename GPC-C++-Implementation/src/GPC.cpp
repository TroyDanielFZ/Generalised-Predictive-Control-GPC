#include "GPC.h"
#include <algorithm>
#include "TMatrix.h"
#include <fstream>


TMatrix GPC::toEplitz(double* dat) {
	_ASSERT(dat != nullptr);
	// Construct the matrix AL/ G from given data
	// given dat = [g0, g1, g2, ... ,g(N-1)]
	// The return is a matrix like
	//
	//     [  g0           0      0      ...   0   ]
	//     [  g1           g0     0      ...   0   ]
	//     [  g2           g1     g0     ...   0   ]
	//     [  .            .      .      ...   .   ]
	//     [  g(Nu)        .      .      ...   g0  ]
	//     [  g(Nu+1)      .      .      ...   g1  ]
	//     [  .            .      .      ...   .   ]
	//     [  g(N-1)       .      .      ...   .   ]

	TMatrix temp(m_N, m_Nu);
	for (int r = 0; r < m_N; ++r) {
		double* data = temp[r];
		int copy_length = std::min(r + 1, m_Nu);
		for(int i = 0; i < copy_length; ++i)
			data[i] = dat[r-i];
	}
	return temp;
}

bool GPC::saveParameters(const std::string& filename) {
	try{
		std::ofstream os(filename, std::ios::binary);
		if (os.bad()) return false;
		os   << "theta  = " << m_Theta  << std::endl
			 << "N      = " << m_N      << std::endl
			 << "Nu     = " << m_Nu     << std::endl
			 << "na     = " << m_na     << std::endl
			 << "nb     = " << m_nb     << std::endl
			 << "rho    = " << m_rho    << std::endl
			 << "alpha  = " << m_alpha  << std::endl
			 << "lambda = " << m_lambda;
	}
	catch (...) {
		return false;
	}
	return true;
}

bool GPC::loadParameters(const std::string& filename) {
	try{
		char ch;
		std::string str;
		std::ifstream in(filename, std::ios::binary);
		if (in.bad()|| !in.is_open()) return false;
		in   >> str >> str >> m_Theta
			 >> str >> str >> m_N
			 >> str >> str >> m_Nu
			 >> str >> str >> m_na
			 >> str >> str >> m_nb
			 >> str >> str >> m_rho
			 >> str >> str >> m_alpha
			 >> str >> str >> m_lambda ;
		// std::cout   << "In GPC:"   << std::endl
		// 		    << "theta  = " << m_Theta   << std::endl
		// 	        << "N      = " << m_N       << std::endl
		// 	        << "Nu     = " << m_Nu      << std::endl
		// 	        << "na     = " << m_na      << std::endl
		// 	        << "nb     = " << m_nb      << std::endl
		// 	        << "rho    = " << m_rho     << std::endl
		// 	        << "alpha  = " << m_alpha   << std::endl
		// 	        << "lambda = " << m_lambda;
	}
	catch (...) {
		return false;
	}
	return true;
}

bool BetaGPC::saveParameters(const std::string& filename) {
	try{
		std::ofstream os(filename, std::ios::binary);
		if (os.bad()) return false;
		os   << "theta  = " << m_Theta  << std::endl
			 << "beta   = " << _beta    << std::endl
			 << "N      = " << m_N      << std::endl
			 << "Nu     = " << m_Nu     << std::endl
			 << "na     = " << m_na     << std::endl
			 << "nb     = " << m_nb     << std::endl
			 << "rho    = " << m_rho    << std::endl
			 << "alpha  = " << m_alpha  << std::endl
			 << "lambda = " << m_lambda;
	}
	catch (...) {
		return false;
	}
	return true;
}

bool BetaGPC::loadParameters(const std::string& filename) {
	try{
		char ch;
		std::string str;
		std::ifstream in(filename, std::ios::binary);
		if (in.bad()|| !in.is_open()) return false;
		in   >> str >> str >> m_Theta
			 >> str >> str >> _beta
			 >> str >> str >> m_N
			 >> str >> str >> m_Nu
			 >> str >> str >> m_na
			 >> str >> str >> m_nb
			 >> str >> str >> m_rho
			 >> str >> str >> m_alpha
			 >> str >> str >> m_lambda ;
		// std::cout  << "In BetaGPC:" << std::endl
		// 		   << "theta  = " << m_Theta   << std::endl
		// 	       << "beta   = " << _beta     << std::endl
		// 	       << "N      = " << m_N       << std::endl
		// 	       << "Nu     = " << m_Nu      << std::endl
		// 	       << "na     = " << m_na      << std::endl
		// 	       << "nb     = " << m_nb      << std::endl
		// 	       << "rho    = " << m_rho     << std::endl
		// 	       << "alpha  = " << m_alpha   << std::endl
		// 	       << "lambda = " << m_lambda;
	}
	catch (...) {
		return false;
	}
	return true;
}
