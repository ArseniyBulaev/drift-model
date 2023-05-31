#pragma once

#include <valarray>

#include "..\Well\Well.h"

namespace MathModel
{
	// ��� �������� ������
	enum class TaskType { BubblesRising };

	class DriftModel
	{
	
	#pragma region Fields
	public:
	const double g = 9.81;  // ��������� ���������� �������
	private:
	TaskType _task_type; // ��������� ������������� ������
	#pragma endregion

	#pragma region Constructor
	public:
	DriftModel(TaskType task_type) : _task_type(task_type) {};
	#pragma endregion

	#pragma region Public Methods
	public:
		// ���������� ����������������� ��������
		double CalculateHydrostaticPressure(double rho, double h);
		// �������� ������ � �������� ������� ����
		double CalculateGasProfileParameter(double alpha_g_0, double U, double d, double p, double dp_dz);
		double CalculateDriftVelocity(double alpha_g_0, double U, double d, double p);
		// �������� �������� ������ � �������� ������� ����
		double CalculateGasProfileParameter_TEST(double alpha_g);
		double CalculateDriftVelocity_TEST(double alpha_g);
		// �������� �����
		double GetMixtureVelocity(double alpha_g, double alpha_l, double v_g, double v_l);
		// �������� ����
		// ���������
		double GetGasDensity(double p);
		double GetLiquidDensity(double p);
		double GetMixtureDensity(double alpha_g, double alpha_l, double p);
		// ��������
		double GetGasViscosity();
		double GetLiquidViscosity();
		double GetMixtureViscosity(double alpha_g);
		// ����������� ������
		double GetFrictionCoefficient(double alpha_g, double alpha_l, double v_m, double p, double d, double eps);
		// ����������� ��������
		// ���������
		double GetCharacteristicGasDensity();
		double GetCharacteristicLiquidDensity();

		// �������� ������ (�������������� ����������)

		// ��������� �������
		double GetGSGasVolumeFractionBoundaryCondition();
		double GetGSPressureBoundaryCondition();
		double GetGSLiquidVelocityBoundaryCondition();
		double GetGSGasVelocityBoundaryCondition();
		double GetGSMixtureVelocityBoundaryCondition();

		// ��������� �������
		double GetGSGasVolumeFractionInitialCondition(double L, double z);


		// �������� ������ (�������������� �������� ��� �������� ��������)

		// ����������� ��������
		double GetPRCharacteristicVelocity(double p_wf, double well_bottom_cross_section_area);
		// ����������� ����� �������� ��������
		double GetPRCharacteristicStopTime(double p_wf, double well_bottom_cross_section_area, double well_length);

		// ��������� �������
		double GetPRLiquidVelocityBoundaryCondition(double p_wf, double well_bottom_cross_section_area);
		double GetPRLiquidVolumeFractionBoundaryCondition(double v_m, double v_l);

		// ��������� �������
		double GetPRGasVolumeFractionInitialCondition();
	#pragma endregion
		
	};
}

