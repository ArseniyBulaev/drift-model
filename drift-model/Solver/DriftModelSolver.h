#pragma once

#include <valarray>
#include <algorithm>

#include "..\Well\Well.h"
#include "..\MathModel\MathModel.h"


class DriftModelSolver
{
public:
	DriftModelSolver(double dz, double dt, const Well& well, MathModel::TaskType task_type);
	void Solve();
	const std::valarray<double> & GetV_m() const;
	const std::valarray<double> & GetV_g() const;
	const std::valarray<double> & GetV_l() const;
	const std::valarray<double> & GetAlpha_g() const;
	const std::valarray<double> & GetAlpha_l() const;
	const std::valarray<double> & Get_P() const;
	const double & GetDt() const;
	const double & GetDz() const;

private:

	const Well _well; // ��������
	const double _dz; // ��� �� ������������
	const double _dt; // ��� �� �������

	MathModel::DriftModel _drift_model; // ����� �������������� ������
	size_t _n_points_cell_properties; // ����� ����� ��� �������� ���� ���������� ����� ���������
	size_t _n_points_cell_velocities; // ����� ����� ��� �������� ���������

	// ������� � ���������� �����������
	std::valarray<double> _theta; // ���� ������� �����
	std::valarray<double> _d; // ������� �����
	std::valarray<double> _eps; // ������������� �������������

	std::valarray<double> _v_m; // ������� �����
	std::valarray<double> _v_g; // ������� ����
	std::valarray<double> _v_l; // ������� ��������
	std::valarray<double> _alpha_g; // �������� ���� ����
	std::valarray<double> _alpha_l; // �������� ���� ��������
	std::valarray<double> _p; // �������� ���������� �����


	void InitializeGeometryParameters(); // ������������� ���������� ��������
	
	int CalculateN();
	
	std::valarray<double> CalculateApproximateMixtureSpeed();
	
	std::valarray<double> CalculateMixtureVelocityCorrection(const std::valarray<double> & p_corr, const std::valarray<double> & v_m_approx);
	
	std::valarray<double> CalculateGasVelocity(const std::valarray<double> & p, const std::valarray<double> & v_m);
	std::valarray<double> CalculateLiquidVelocity(const std::valarray<double> & alpha_g, const std::valarray<double> & alpha_l, const std::valarray<double> & v_g, const std::valarray<double> & v_m);
	
	std::valarray<double> CalculateGasVolumeFraction(const std::valarray<double>& p_current, const std::valarray<double>& v_g);
	std::valarray<double> CalculateLiquidVolumeFraction(const std::valarray<double>& p_current, const std::valarray<double>& v_l);
	std::valarray<double> CalculatePressureCorrection(const std::valarray<double> & v_m_approx, const std::valarray<double> & alpha_g_past, const std::valarray<double> & alpha_l_past, const std::valarray<double> & p_past);

	// TEST
	std::valarray<double> CalculateGasVelocity_TEST(const std::valarray<double> & v_m);

	double CalculateGasImbalance(const std::valarray<double>& alpha_g, const std::valarray<double>& alpha_g_past, const std::valarray<double>& p_g, const std::valarray<double>& p_g_past, const std::valarray<double>& v_g);
	double CalculateLiquidImbalance(const std::valarray<double>& alpha_l, const std::valarray<double>& alpha_l_past, const std::valarray<double>& p_l, const std::valarray<double>& p_l_past, const std::valarray<double>& v_l);

	void TDMA(std::valarray<double> & v, const std::valarray<double> & a, const std::valarray<double> & b, const std::valarray<double> & c, const std::valarray<double> & d);
	void SIMPLE();
};

