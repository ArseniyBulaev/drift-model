#pragma once

#include <valarray>
#include <algorithm>

#include "..\Well\Well.h"
#include "..\MathModel\MathModel.h"
#include"..\Writer.h"


class DriftModelSolver
{
public:
	DriftModelSolver(double dz, double dt, const Well& well, MathModel::TaskType task_type);
	void Solve();
	const std::valarray<double> & GetV_m() const;
	const std::valarray<double> & GetV_g() const;
	const std::valarray<double> & GetV_l() const;
	const std::valarray<double> & GetAlpha_g() const;
	const std::valarray<double> & Get_P() const;
	const double & GetDt() const;
	const double & GetDz() const;

private:

	const Well _well; // ��������
	const double _dz; // ��� �� ������������
	const double _dt; // ��� �� �������

	MathModel::DriftModel _drift_model; // ����� �������������� ������
	Writer _results_writer; // ������ ��� ������ �����������

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
	std::valarray<double> _p; // �������� ���������� �����


	void InitializeGeometryParameters(); // ������������� ���������� ��������
	
	int CalculateNumberOfPoints();
	
	void CalculateApproximateMixtureSpeed(std::valarray<double> & v_m_intermediate);
	
	std::valarray<double> CalculateMixtureVelocityCorrection(const std::valarray<double> & p_corr);
	
	std::valarray<double> CalculateGasVelocity(const std::valarray<double> & p, const std::valarray<double> & v_m);
	
	std::valarray<double> CalculateGasVolumeFraction(const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g);
	std::valarray<double> CalculatePressureCorrection(const std::valarray<double> & v_m_intermediate);

	// TEST
	std::valarray<double> CalculateGasVelocity_TEST(const std::valarray<double> & v_m);

	double CalculateGasImbalance(const std::valarray<double>& alpha_g_intermediate, const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g_intermediate);
	

	void TDMA(std::valarray<double> & v, const std::valarray<double> & a, const std::valarray<double> & b, const std::valarray<double> & c, const std::valarray<double> & d);
	void SimpleAlgorithm();
};

