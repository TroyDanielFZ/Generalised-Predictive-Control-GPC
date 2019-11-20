#pragma once
#include <cmath>
#include "TMath.h"
#include "TMatrix.h"
#include <iostream>
#include <string>
class GPC {
protected:
	int m_N, m_Nu, m_na, m_nb;
	double m_rho, m_alpha, m_lambda;
	TCircularBuffer  m_Y, m_U, m_Theta, m_A, m_B;
	TMatrix m_P;

	void Diophantine(TMatrix& F, TMatrix& H, TCircularBuffer& g,
		const TCircularBuffer& A, const TCircularBuffer& B, int N) {
		double* e = new double[N] {0};

		e[0] = 1;
		g.push(B[0]);
		F.Row(0) = -A.slice(1, -1); // truncate the first item
		H.Row(0) = B.slice(1, -1); // truncate the first item
		auto F0 = F.Row(0).toCircularBuffer();
		auto H0 = H.Row(0).toCircularBuffer();
		for (int i = 1; i < N; ++i) {
			e[i] = F[i - 1][0];
			F.Row(i) = F.Row(i - 1).toCircularBuffer().push(0) + (F0 * e[i]);

			g.push(e[i] * B[0] + H[i - 1][0]);
			H.Row(i) = H.Row(i - 1).toCircularBuffer().push(0) + H0 * e[i];
		}
		delete[] e;
	}
public:
	TMatrix toEplitz(double* dat);
	GPC(int na, int nb, int N, int Nu, double rho, double alpha, double lambda)
		:m_na(na)
		, m_nb(nb)
		, m_N(N)
		, m_Nu(Nu)
		, m_rho(rho)
		, m_alpha(alpha)
		, m_lambda(lambda)
		, m_Y(na + 1)
		, m_U(nb + 2)
		, m_P(TMatrix::createIdentity(na + nb + 1) * 100)
		, m_Theta(na + nb + 1)
		, m_A(na + 1)
		, m_B(nb + 1) {
		//m_pA = new double[na + 1L]{ 0.0f };
		//m_pB = new double[nb + 2L]{ 0.0f };
		double init = 1.0;
		for (int i = 0; i < na + nb + 1; ++i)
			m_Theta.push((init *= 0.7));
	}
	~GPC() { }

	virtual double getControl(double yr, double y) {
		TMatrix Phi = TMatrix::augment(m_Y.toReverse().differ(), m_U.toReverse().differ(false), false);
		double deltaY = y - m_Y.last();
		m_Y.push(y);
		auto PhiT = Phi.transpose();
		double ratio = (Phi.transpose() * m_P * Phi)[0][0] + m_rho;
		TMatrix mat = m_P * Phi * (deltaY - PhiT * m_Theta.toMatrix());
		mat /= ratio;

		// Update the estimations
		m_Theta += TCircularBuffer(mat); // I'm pretty sure that mat is a column vector
		m_P = (m_P - m_P * Phi * PhiT * m_P / ratio) / m_rho;

		m_A.push(1);
		for (int i = 0; i < m_na; ++i) m_A.push(m_Theta[i]);
		for (int i = m_na; i <= m_na + m_nb; ++i) m_B.push(m_Theta[i]);

		TCircularBuffer g(m_N);
		TMatrix F(m_N, m_na + 1), H(m_N, m_nb);
		TCircularBuffer A(m_na + 2);
		A.push(1);
		for (int i = 0; i < m_na; ++i)
			A.push(m_A[i + 1] - m_A[i]);
		A.push(-m_A[m_na]);

		Diophantine(F, H, g, A, m_B, m_N);

		auto Y0 = F * m_Y.toReverse().toMatrix() - H * m_U.slice(1, -1).toReverse().differ();
		double * ptr_g = new double[g.length()]{ 0.0 };
		for (long i = 0; i < g.length(); ++i) {
			ptr_g[i] = g[i];
		}
		TMatrix G = toEplitz(ptr_g);
		delete[] ptr_g;

		auto GT = G.transpose();
		TMatrix gain = ( GT * G + m_lambda * TMatrix::createIdentity(m_Nu)).divideBy(GT).slice(0, 0, 0,-1);
		TCircularBuffer Yd(m_N);
		// generate the tracking trajectory
		double dY = y - yr;
		for (int i = 0; i < m_N; ++i) {
			dY *= m_alpha;
			Yd.push(yr + dY);
		}
		double deltaU =  (gain * (Yd.toMatrix() - Y0))[0][0];

		m_U.push(m_U.last() + deltaU);

		return m_U.last();

	}
	virtual bool saveParameters(const std::string & filename);
	virtual bool loadParameters(const std::string& filename);
};

class BetaGPC : public GPC {
private:
	double _beta;
public:
	BetaGPC(int na, int nb, int N, int Nu, double rho, double alpha, double lambda, double beta)
		:GPC(na, nb, N, Nu, rho, alpha, lambda)
		, _beta(beta) { }


	double getControl(double yr, double y) {
		using namespace std;
		TMatrix Phi = TMatrix::augment(m_Y.toReverse().differ(), m_U.toReverse().differ(false), false);
		double deltaY = y - m_Y.last();
		m_Y.push(y);
		auto PhiT = Phi.transpose();
		double ratio = (Phi.transpose() * m_P * Phi)[0][0] + m_rho;
		TMatrix mat = m_P * Phi * (deltaY - PhiT * m_Theta.toMatrix());
		mat /= ratio;

		// Update the estimations
		m_Theta += TCircularBuffer(mat); // I'm pretty sure that mat is a column vector
		m_P = (m_P - m_P * Phi * PhiT * m_P / ratio) / m_rho;

		m_A.push(1);
		for (int i = 0; i < m_na; ++i) m_A.push(m_Theta[i]);
		for (int i = m_na; i <= m_na + m_nb; ++i) m_B.push(m_Theta[i]);

		TCircularBuffer g(m_N);
		TMatrix F(m_N, m_na + 1), H(m_N, m_nb);
		TCircularBuffer A(m_na + 2);
		// A = conv(m_A, [1, -1]);
		A.push(1);
		for (int i = 0; i < m_na; ++i)
			A.push(m_A[i + 1] - m_A[i]);
		A.push(-m_A[m_na]);

		Diophantine(F, H, g, A, m_B, m_N);

		auto Y0 = F * m_Y.toReverse().toMatrix() - H * m_U.slice(1, -1).toReverse().differ();
		double * ptr_g = new double[g.length()]{ 0.0 };
		for (long i = 0; i < g.length(); ++i) {
			ptr_g[i] = g[i];
		}
		TMatrix G = toEplitz(ptr_g);
		delete[] ptr_g;
		auto GT = G.transpose();
		TMatrix gain = ( GT * G + m_lambda * TMatrix::createIdentity(m_Nu)).divideBy(GT).slice(0, 0, 0,-1);
		TCircularBuffer Yd(m_N);
		// generate the tracking trajectory
		double dY = y - yr;
		for (int i = 0; i < m_N; ++i) {
			dY *= m_alpha;
			Yd.push(yr + dY);
		}
		double deltaU = _beta * (gain * (Yd.toMatrix() - Y0))[0][0];

		m_U.push(m_U.last() + deltaU);

		return m_U.last();
	}
	bool saveParameters(const std::string& filename);
	bool loadParameters(const std::string& filename);
};
