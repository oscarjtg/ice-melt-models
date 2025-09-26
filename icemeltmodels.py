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
        # The seawater freezing temperature at the given temperature, salinity, and pressure.
        freezing_temperature = self.liquidus_slope * salinity + self.liquidus_intercept + self.liquidus_pressure_coefficient * pressure

        # The constant of proportionality.
        constant = (self.density_water / self.density_ice) * (self.heat_capacity_water / self.latent_heat_ice) * self.transfer_coefficient

        # Return the melt rate.
        return constant * current_speed * (temperature - freezing_temperature)
    

class TwoEquationMeltModel(MeltModel):
    def __init__(self, transfer_coefficient=5.909314681e-04, ice_temperature=-10, **kwargs):
        super().__init__(**kwargs)
        self.transfer_coefficient = transfer_coefficient
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
        # The seawater freezing temperature at the given temperature, salinity, and pressure.
        freezing_temperature = self.liquidus_slope * salinity + self.liquidus_intercept + self.liquidus_pressure_coefficient * pressure

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


class ThreeEquationMeltModel(MeltModel):
    def __init__(self, heat_transfer_coefficient=1.083374358e-03, salt_transfer_coefficient=3.053145919e-05, ice_temperature=-10, **kwargs):
        super().__init__(**kwargs)
        self.heat_transfer_coefficient = heat_transfer_coefficient
        self.salt_transfer_coefficient = salt_transfer_coefficient
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
        # Solve quadratic equation for the melt rate.
        alpha = (self.density_water / self.density_ice) * self.salt_transfer_coefficient * current_speed
        beta = self.density_water * self.heat_capacity_water * self.heat_transfer_coefficient * current_speed
        gamma = self.density_ice * (self.latent_heat_ice - self.heat_capacity_ice * self.ice_temperature)
        delta = self.density_ice * self.heat_capacity_ice
        epsilon = self.liquidus_intercept + self.liquidus_pressure_coefficient * pressure
        
        # Coefficients of the quadratic (a, b, and c).
        a = gamma + delta * epsilon
        b = alpha * gamma - beta * temperature + alpha * delta * self.liquidus_slope * salinity + alpha * delta * epsilon + beta * epsilon
        c = -alpha * beta * temperature + alpha * beta * self.liquidus_slope * salinity + alpha * beta * epsilon

        # Quadratic formula.
        melt_rate = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
        return melt_rate
        