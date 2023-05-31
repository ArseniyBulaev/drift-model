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
		double CalculateHydrostaticPressure(double rho, double h);
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
		// Характерные величины
		// Плотность
		double GetCharacteristicGasDensity();
		double GetCharacteristicLiquidDensity();

		// Тестовая задача (Гравитационная сегрегация)

		// Граничные условия
		double GetGSGasVolumeFractionBoundaryCondition();
		double GetGSPressureBoundaryCondition();
		double GetGSLiquidVelocityBoundaryCondition();
		double GetGSGasVelocityBoundaryCondition();
		double GetGSMixtureVelocityBoundaryCondition();

		// Начальные условия
		double GetGSGasVolumeFractionInitialCondition(double L, double z);


		// Тестовая задача (Восстановление давления при закрытии скважины)

		// Характерная скорость
		double GetPRCharacteristicVelocity(double p_wf, double well_bottom_cross_section_area);
		// Характерное время закрытия скважины
		double GetPRCharacteristicStopTime(double p_wf, double well_bottom_cross_section_area, double well_length);

		// Граничные условия
		double GetPRLiquidVelocityBoundaryCondition(double p_wf, double well_bottom_cross_section_area);
		double GetPRLiquidVolumeFractionBoundaryCondition(double v_m, double v_l);

		// Начальные условия
		double GetPRGasVolumeFractionInitialCondition();
	#pragma endregion
		
	};
}

