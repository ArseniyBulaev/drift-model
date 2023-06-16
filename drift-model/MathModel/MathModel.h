#pragma once

#include <valarray>

#include "..\Well\Well.h"

namespace MathModel
{
	// Тип Решаемой задачи
	enum class TaskType { BubblesRising, Debug };

	class DriftModel
	{
	
		#pragma region Fields
		public:
		const double g = 9.81;  // Ускорение свободного падения
		const double atm = 101300; // Атмосферное давление
		private:
		TaskType _task_type; // Выбранная пользователем задача
		double _U = 0; // Характерная скорость
		double _L = 0; // Характерная длина
		double _alpha_g_0 = 0; // Характерная объёмная доля пузырьков
		#pragma endregion

		#pragma region Constructor
		public:
		DriftModel(TaskType task_type, const Well & well) : _task_type(task_type), _L(well.GetLength())
		{
			switch (task_type)
			{
			case MathModel::TaskType::BubblesRising:
			{
				_alpha_g_0 = 0.05;
				SetBubblesRisingCharacteristicVelocity(well);
				break;
			}
			case MathModel::TaskType::Debug:
				SetDebugCharacteristicVelocity(well);
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
			const std::valarray<double>& theta,
			double dz);
		void SetBoundaryConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			const Well well,
			double dt);

	
		// Скорость дрейфа и параметр профиля газа
		std::valarray<double> CalculateC_0(const std::valarray<double> & d, const std::valarray<double> & alpha_g,  const std::valarray<double>& p, size_t n_points_cell_velocities, double dz);
		std::valarray<double> CalculateV_d(const std::valarray<double> & d, const std::valarray<double> & alpha_g, std::valarray<double> p, size_t n_points_cell_velocities);

		double GasSteadyFlowAnalyticsVelocity(double p, double z, double F);
		double GasSteadyFlowAnalyticsDensity(double v, double p);

		// Скорость смеси
		double GetMixtureVelocity(double alpha_g, double alpha_l, double v_g, double v_l);
		// Скорость воды
		// Плотность
		double GetGasDensity(double p, double z);
		double GetLiquidDensity(double p);
		double GetMixtureDensity(double alpha_g, double alpha_l, double p, double z);
		// Вязкость
		double GetLiquidViscosity();
		double GetMixtureViscosity(double alpha_g);
		// Коэффициент трения
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
			const std::valarray<double>& theta,
			double dz);

		void SetBubblesRisingBoundaryConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			const Well well,
			double dt);

		void SetDebugInitialConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			double dz);

		void SetDebugBoundaryConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			const Well well,
			double dt);

		double GetBubblesRisingLiquidFlow();
		double GetBubblesRisingGasFlow();
		void SetBubblesRisingCharacteristicVelocity(Well well);
		void SetDebugCharacteristicVelocity(Well well);
		void SetDebugCharacteristicGasVolumeFraction();
		#pragma endregion
		

		



		double CalculateHydrostaticPressure(double rho, double h);
		#pragma endregion

		
	};
}

