###########################################################
#                                                         #
# icemeltmodels.py                                        #
#                                                         #
# This module contains classes that provide different     #
# models for melt at an ice-ocean interface.              #
#                                                         #
# Author: Oscar Tovey Garcia                              #
#                                                         #
###########################################################

import numpy as np

class MeltModel():
    def __init__(self, 
                 density_ice=916.0, 
                 heat_capacity_ice=2000.0,
                 latent_heat_ice=3.34e05, 
                 density_water=1030.0,
                 heat_capacity_water=3974.0, 
                 liquidus_slope=-0.0573, 
                 liquidus_intercept=0.0832, 
                 liquidus_pressure_coefficient=-7.53e-08
                 ):
        self.density_ice = density_ice
        self.heat_capacity_ice = heat_capacity_ice
        self.latent_heat_ice = latent_heat_ice
        self.density_water = density_water
        self.heat_capacity_water = heat_capacity_water
        self.liquidus_slope = liquidus_slope
        self.liquidus_intercept = liquidus_intercept
        self.liquidus_pressure_coefficient = liquidus_pressure_coefficient

    def freezing_temperature(self, current_speed, temperature, salinity, pressure):
        """
        Calculates the seawater freezing point temperature according to the 
        linear equation of state from Jenkins et al (2010).

        Parameters:
            current_speed: The magnitude of the velocity of the free-stream current
                           adjacent to the ice interface but outside the turbulent 
                           boundary layer, in m/s. (Not used)

            temperature:   The temperature of the ocean water, in degrees Celsius. (Not used)

            salinity:      The salinity of the ocean water, in g/kg.

            pressure:      The pressure at that point, in Pa.
        
        Returns:
            The freezing temperature predicted by this model.
        """
        # The seawater freezing temperature at the given ambient salinity and pressure.
        freezing_temperature = self.liquidus_slope * salinity + self.liquidus_intercept + self.liquidus_pressure_coefficient * pressure
        return freezing_temperature


class TwoEquationMeltModelNeglectingConduction(MeltModel):
    def __init__(self, transfer_coefficient=5.909314681e-04, **kwargs):
        super().__init__(**kwargs)
        self.transfer_coefficient = transfer_coefficient

    def melt_rate(self, current_speed, temperature, salinity, pressure):
        """
        Calculates the melt rate according to the two equation melt model 
        from Jenkins et al (2010).

        Turbulent heat flux from the ocean to the ice is parameterised 
        by assuming that this flux is proportional to both the speed of 
        the free-stream current and also the difference between 
        the ocean temperature and the seawater freezing temperature.

        Conductive heat flux into the ice is assumed to be negligible.

        Parameters:
            current_speed: The magnitude of the velocity of the free-stream current
                           adjacent to the ice interface but outside the turbulent 
                           boundary layer, in m/s.

            temperature:   The temperature of the ocean water, in degrees Celsius.

            salinity:      The salinity of the ocean water, in g/kg.

            pressure:      The pressure at that point, in Pa.

        Returns:
            The melt rate predicted by this model.
        """
        # The seawater freezing temperature at the given ambient salinity and pressure.
        freezing_temperature = self.freezing_temperature(current_speed, temperature, salinity, pressure)

        # The constant of proportionality.
        constant = (self.density_water / self.density_ice) * (self.heat_capacity_water / self.latent_heat_ice) * self.transfer_coefficient

        # Return the melt rate.
        return constant * current_speed * (temperature - freezing_temperature)
    

class TwoEquationMeltModel(TwoEquationMeltModelNeglectingConduction):
    def __init__(self, ice_temperature=-10, **kwargs):
        super().__init__(**kwargs)
        self.ice_temperature = ice_temperature

    def melt_rate(self, current_speed, temperature, salinity, pressure):
        """
        Calculates the melt rate according to the two equation melt model 
        from Jenkins et al (2010).

        Turbulent heat flux from the ocean to the ice is parameterised 
        by assuming that this flux is proportional to both the speed of 
        the free-stream current and also the difference between 
        the ocean temperature and the seawater freezing temperature.

        Conductive heat flux into the ice is parameterised by assuming that this 
        flux is proportional to both the melt rate and also to the difference 
        between the ice temperature and the seawater freezing temperature.

        Parameters:
            current_speed: The magnitude of the velocity of the free-stream current
                           adjacent to the ice interface but outside the turbulent 
                           boundary layer, in m/s.

            temperature:   The temperature of the ocean water, in degrees Celsius.

            salinity:      The salinity of the ocean water, in g/kg.

            pressure:      The pressure at that point, in Pa.

        Returns:
            The melt rate predicted by this model.
        """
        # The seawater freezing temperature at the given ambient salinity and pressure.
        freezing_temperature = self.freezing_temperature(current_speed, temperature, salinity, pressure)

        # The constant of proportionality.
        numerator = self.density_water * self.heat_capacity_water * self.transfer_coefficient
        denominator = self.density_ice * (self.latent_heat_ice + self.heat_capacity_ice * (freezing_temperature - self.ice_temperature))
        constant = numerator / denominator

        # Return the melt rate.
        return constant * current_speed * (temperature - freezing_temperature)
    

class ThreeEquationMeltModelNeglectingConduction(MeltModel):
    def __init__(self, heat_transfer_coefficient=1.083374358e-03, salt_transfer_coefficient=3.053145919e-05, **kwargs):
        super().__init__(**kwargs)
        self.heat_transfer_coefficient = heat_transfer_coefficient
        self.salt_transfer_coefficient = salt_transfer_coefficient

    def melt_rate(self, current_speed, temperature, salinity, pressure):
        """
        Calculates the melt rate according to the three equation melt model 
        from Jenkins et al (2010).

        Turbulent heat flux from the ocean to the ice is parameterised 
        by assuming that this flux is proportional to both the speed of 
        the free-stream current and also the difference between 
        the ocean temperature and the seawater freezing temperature.

        Conductive heat flux into the ice is assumed to be negligible.

        Parameters:
            current_speed: The magnitude of the velocity of the free-stream current
                           adjacent to the ice interface but outside the turbulent 
                           boundary layer, in m/s.

            temperature:   The temperature of the ocean water, in degrees Celsius.

            salinity:      The salinity of the ocean water, in g/kg.

            pressure:      The pressure at that point, in Pa.

        Returns:
            The melt rate predicted by this model.
        """
        # Solve quadratic equation for the melt rate.
        # a, b, and c are the coefficients of the quadratic equation.
        a = 1.0
        b = (self.density_water / self.density_ice) * self.salt_transfer_coefficient * current_speed + (self.density_water / self.density_ice) * (self.heat_capacity_water / self.latent_heat_ice) * self.heat_transfer_coefficient * current_speed * (self.liquidus_intercept + self.liquidus_pressure_coefficient * pressure - temperature)
        c = (self.density_water / self.density_ice)**2 * (self.heat_capacity_water / self.latent_heat_ice) * self.heat_transfer_coefficient * self.salt_transfer_coefficient * current_speed**2 * (self.liquidus_slope * salinity + self.liquidus_intercept + self.liquidus_pressure_coefficient * pressure - temperature)

        # Quadratic formula.
        melt_rate = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
        return melt_rate
    
    def boundary_salinity(self, current_speed, temperature, salinity, pressure):
        """
        Calculates the salinity at the ice-ocean boundary according to the 
        three equation melt model from Jenkins et al (2010), 
        neglecting conductive heat flux into the ice.

        Parameters:
            current_speed: The magnitude of the velocity of the free-stream current
                           adjacent to the ice interface but outside the turbulent 
                           boundary layer, in m/s.

            temperature:   The temperature of the ocean water, in degrees Celsius.

            salinity:      The salinity of the ocean water, in g/kg.

            pressure:      The pressure at that point, in Pa.

        Returns:
            The salinity at the ice-ocean boundary.
        """
        melt_rate = self.melt_rate(current_speed, temperature, salinity, pressure)
        numerator = self.density_water * self.salt_transfer_coefficient * current_speed
        denominator = self.density_ice * melt_rate + numerator
        constant = np.where(denominator == 0.0, np.nan, numerator / denominator)

        return constant * salinity
    
    def boundary_temperature(self, current_speed, temperature, salinity, pressure):
        """
        Calculates the temperature at the ice-ocean boundary according to the 
        three equation melt model from Jenkins et al (2010), 
        neglecting conductive heat flux into the ice.

        Parameters:
            current_speed: The magnitude of the velocity of the free-stream current
                           adjacent to the ice interface but outside the turbulent 
                           boundary layer, in m/s.

            temperature:   The temperature of the ocean water, in degrees Celsius.

            salinity:      The salinity of the ocean water, in g/kg.

            pressure:      The pressure at that point, in Pa.

        Returns:
            The temperature at the ice-ocean boundary.
        """
        boundary_salinity = self.boundary_salinity(current_speed, temperature, salinity, pressure)

        return self.freezing_temperature(current_speed, temperature, boundary_salinity, pressure)

class ThreeEquationMeltModel(ThreeEquationMeltModelNeglectingConduction):
    def __init__(self, ice_temperature=-10, **kwargs):
        super().__init__(**kwargs)
        self.ice_temperature = ice_temperature
    
    def melt_rate(self, current_speed, temperature, salinity, pressure):
        """
        Calculates the melt rate according to the three equation melt model 
        from Jenkins et al (2010).

        Turbulent heat flux from the ocean to the ice is parameterised 
        by assuming that this flux is proportional to both the speed of 
        the free-stream current and also the difference between 
        the ocean temperature and the seawater freezing temperature.

        Conductive heat flux into the ice is parameterised by assuming that this 
        flux is proportional to both the melt rate and also to the difference 
        between the ice temperature and the seawater freezing temperature.

        Parameters:
            current_speed: The magnitude of the velocity of the free-stream current
                           adjacent to the ice interface but outside the turbulent 
                           boundary layer, in m/s.

            temperature:   The temperature of the ocean water, in degrees Celsius.

            salinity:      The salinity of the ocean water, in g/kg.

            pressure:      The pressure at that point, in Pa.

        Returns:
            The melt rate predicted by this model.
        """
        rho_w = self.density_water
        rho_i = self.density_ice
        c_w = self.heat_capacity_water
        c_i = self.heat_capacity_ice
        L_i = self.latent_heat_ice
        T_i = self.ice_temperature

        gamma_T = self.heat_transfer_coefficient
        gamma_S = self.salt_transfer_coefficient
        
        lambda_1 = self.liquidus_slope
        lambda_2 = self.liquidus_intercept
        lambda_3 = self.liquidus_pressure_coefficient

        U = current_speed # Not used here, since S_b is independent of U.
        T_w = temperature
        S_w = salinity
        P_b = pressure
        
        # Coefficients of the quadratic (a, b, and c).
        a = (L_i - c_i * (T_i - lambda_2 - lambda_3 * P_b)) * rho_i**2
        b = (gamma_S * L_i 
             + gamma_S * c_i * (lambda_1 * S_w + lambda_2 + lambda_3 * P_b - T_i)
             - gamma_T * c_w * (T_w - lambda_2 - lambda_3 * P_b)
             ) * rho_w * rho_i * U
        c = (lambda_1 * S_w + lambda_2 + lambda_3 * P_b - T_w) * gamma_T * gamma_S * rho_w**2 * c_w * U**2

        # Quadratic formula.
        melt_rate = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
        return melt_rate
    
    def boundary_salinity(self, current_speed, temperature, salinity, pressure):
        """
        Solves a quadratic equation for the boundary salinity found by eliminating 
        the melt rate and the boundary temperature from the Three Equation melt model
        in Jenkins et al. (2010)
        
        Parameters:
            current_speed: The magnitude of the velocity of the free-stream current
                           adjacent to the ice interface but outside the turbulent 
                           boundary layer, in m/s.

            temperature:   The temperature of the ocean water, in degrees Celsius.

            salinity:      The salinity of the ocean water, in g/kg.

            pressure:      The pressure at that point, in Pa.

        Returns:
            The boundary salinity predicted by this model.
        """
        c_w = self.heat_capacity_water
        c_i = self.heat_capacity_ice
        L_i = self.latent_heat_ice
        T_i = self.ice_temperature

        gamma_T = self.heat_transfer_coefficient
        gamma_S = self.salt_transfer_coefficient
        
        lambda_1 = self.liquidus_slope
        lambda_2 = self.liquidus_intercept
        lambda_3 = self.liquidus_pressure_coefficient

        U = current_speed # Not used here, since S_b is independent of U.
        T_w = temperature
        S_w = salinity
        P_b = pressure

        # Coefficients of quadratic equation.
        a = (gamma_S * c_i - gamma_T * c_w) * lambda_1  # in m^2 s^-2 (g/kg)^-1
        b = (gamma_S * L_i 
             + gamma_T * c_w * (T_w - lambda_2 - lambda_3 * P_b)
             - gamma_S * c_i * (lambda_1 * S_w + T_i - lambda_2 - lambda_3 * P_b) # in m^2 s^-2
        )
        c = (c_i * (T_i - lambda_2 - lambda_3 * P_b) - L_i) * gamma_S * S_w # in m^2 s^-2 g/kg

        # Solve quadratic equation. Return the positive root.
        boundary_salinity = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a) # in g/kg
        return boundary_salinity * (U / U) # * (U / U) is to have array shape of U.
        