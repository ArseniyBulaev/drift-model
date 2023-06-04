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
		double _U = 0; // ����������� ��������
		double _L = 0; // ����������� �����
		double _alpha_g_0 = 0; // ����������� �������� ���� ���������
		#pragma endregion

		#pragma region Constructor
		public:
		DriftModel(TaskType task_type, const Well & well) : _task_type(task_type), _L(well.GetLength())
		{
			switch (task_type)
			{
			case MathModel::TaskType::BubblesRising:
				SetBubblesRisingCharacteristicVelocity(well);
				SetBubblesRisingCharacteristicGasVolumeFraction();
				break;
			}
		};
		#pragma endregion

		#pragma region Public Methods
		public:
		void SetInitialConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			double dz);
		void SetBoundaryConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			const Well well,
			double dt);

	
		// �������� ������ � �������� ������� ����
		std::valarray<double> CalculateC_0(const std::valarray<double> & d, const std::valarray<double>& p, size_t n_points_cell_velocities, double dz);
		std::valarray<double> CalculateV_d(const std::valarray<double> & d, std::valarray<double> p, size_t n_points_cell_velocities);


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
		double GetLiquidViscosity();
		double GetMixtureViscosity(double alpha_g);
		// ����������� ������
		double GetFrictionCoefficient(double alpha_g, double alpha_l, double v_m, double p, double d, double eps);
	
		#pragma region Characteristic Values
		double GetCharacteristicGasDensity();
		double GetCharacteristicLiquidDensity();
		#pragma endregion

	
		#pragma endregion

		#pragma region Private Methods
		private:
		#pragma region Bubbles Rising Task

		void SetBubblesRisingInitialConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			double dz);

		void SetBubblesRisingBoundaryConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			const Well well,
			double dt);

		double GetBubblesRisingLiquidFlow(double dt);
		double GetBubblesRisingGasFlow(double dt);
		void SetBubblesRisingCharacteristicVelocity(Well well);
		void SetBubblesRisingCharacteristicGasVolumeFraction();
		#pragma endregion
		

		double CalculateHydrostaticPressure(double rho, double h);
		#pragma endregion

		
	};
}

