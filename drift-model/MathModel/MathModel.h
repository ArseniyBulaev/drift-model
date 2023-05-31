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
		// Вычисление гидростатического давления
		void SetInitialConditions(std::valarray<double>& alpha_g, std::valarray<double>& p, std::valarray<double>& v_m, double dz);

	
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

