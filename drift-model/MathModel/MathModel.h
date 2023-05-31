#pragma once

#include <valarray>

#include "..\Well\Well.h"

namespace MathModel
{
	// Тип Решаемой задачи
	enum class TaskType { BubblesRising };

	class DriftModel
	{
	
		#pragma region Fields
		public:
		const double g = 9.81;  // Ускорение свободного падения
		private:
		TaskType _task_type; // Выбранная пользователем задача
		#pragma endregion

		#pragma region Constructor
		public:
		DriftModel(TaskType task_type) : _task_type(task_type) {};
		#pragma endregion

		#pragma region Public Methods
		public:
		void SetInitialConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& alpha_l,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			double dz);
		void SetBoundaryConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& alpha_l,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			const Well well,
			double t);

	
		// Скорость дрейфа и параметр профиля газа
		double CalculateGasProfileParameter(double alpha_g_0, double U, double d, double p, double dp_dz);
		double CalculateDriftVelocity(double alpha_g_0, double U, double d, double p);
		// Тестовая скорость дрейфа и параметр профиля газа
		double CalculateGasProfileParameter_TEST(double alpha_g);
		double CalculateDriftVelocity_TEST(double alpha_g);
		// Скорость смеси
		double GetMixtureVelocity(double alpha_g, double alpha_l, double v_g, double v_l);
		// Скорость воды
		// Плотность
		double GetGasDensity(double p);
		double GetLiquidDensity(double p);
		double GetMixtureDensity(double alpha_g, double alpha_l, double p);
		// Вязкость
		double GetGasViscosity();
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
			std::valarray<double>& alpha_l,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			double dz);
		void SetBubblesRisingBoundaryConditions(
			std::valarray<double>& alpha_g,
			std::valarray<double>& alpha_l,
			std::valarray<double>& p,
			std::valarray<double>& v_m,
			std::valarray<double>& v_g,
			std::valarray<double>& v_l,
			const Well well,
			double t);
		double GetBubblesRisingLiquidFlow(double t);
		double GetBubblesRisingGasFlow(double t);
		#pragma endregion

		

		double CalculateHydrostaticPressure(double rho, double h);
		#pragma endregion

		
	};
}

