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
		void SetInitialConditions(std::valarray<double>& alpha_g, std::valarray<double>& p, std::valarray<double>& v_m, double dz);

	
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
	
		#pragma region CharacteristicValues
		double GetCharacteristicGasDensity();
		double GetCharacteristicLiquidDensity();
		#pragma endregion

	
		#pragma endregion

		#pragma region Private Methods
		private:
		void SetBubblesRisingInitialConditions(std::valarray<double>& alpha_g, std::valarray<double>& p, std::valarray<double>& v_m, double dz);
		double CalculateHydrostaticPressure(double rho, double h);
		#pragma endregion

		
	};
}

