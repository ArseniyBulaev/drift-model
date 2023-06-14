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

private:

	const Well _well; // ��������
	const double _dz; // ��� �� ������������
	double _dt; // ��� �� �������

	MathModel::DriftModel _drift_model; // ����� �������������� ������
	Writer _results_writer; // ������ ��� ������ �����������

	size_t _n_points_cell_properties; // ����� ����� ��� �������� ���� ���������� ����� ���������
	size_t _n_points_cell_velocities; // ����� ����� ��� �������� ���������

	double alpha_p_relax = 0.9; // ����������� ���������� ��� ��������
	double alpha_v_relax = 0.3; // ����������� ���������� ��� ��������

	// ������� � ���������� �����������
	std::valarray<double> _theta; // ���� ������� �����
	std::valarray<double> _d; // ������� �����
	std::valarray<double> _eps; // ������������� �������������

	std::valarray<double> _v_m; // ������� �����
	std::valarray<double> _v_g; // ������� ����
	std::valarray<double> _v_l; // ������� ��������
	std::valarray<double> _alpha_g; // �������� ���� ����
	std::valarray<double> _p; // �������� ���������� �����


	void CorrectTimeStep(const std::valarray<double>& v_m_intermediate);
	void InitializeGeometryParameters(); // ������������� ���������� ��������
	int CalculateNumberOfPoints();
	void CalculateApproximateMixtureVelocity(std::valarray<double> & v_m_intermediate);
	std::valarray<double> CalculateMixtureVelocityCorrection(const std::valarray<double> & p_corr);
	std::valarray<double> CalculateGasVelocity(const std::valarray<double> & p, const std::valarray<double> & v_m);
	std::valarray<double> CalculateLiquidVelocity(const std::valarray<double>& v_m, const std::valarray<double>& alpha_g, const std::valarray<double>& v_g);
	std::valarray<double> CalculateGasVolumeFraction(const std::valarray<double>& p_intermediate, const std::valarray<double>& v_g);
	std::valarray<double> CalculatePressureCorrection(const std::valarray<double> & v_m_intermediate, const std::valarray<double>& p_intermediate);
	std::valarray<double> CalculateGasVelocityGeneral(const std::valarray<double> & v_m);
	double CalculateGasImbalance(const std::valarray<double>& alpha_g_intermediate,
		const std::valarray<double>& p_intermediate,
		const std::valarray<double>& v_g_intermediate,
		const std::valarray<double>& alpha_g_previous_iteration,
		const std::valarray<double>& p_previous_iteration,
		const std::valarray<double>& v_g_previous_iteration);
	void TDMA(std::valarray<double> & v, const std::valarray<double> & a, const std::valarray<double> & b, const std::valarray<double> & c, const std::valarray<double> & d);
	void SimpleAlgorithm();
};

