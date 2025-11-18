#
# UtahLSM
#
# Copyright (c) 2017–2025 Jeremy A. Gibbs
# Copyright (c) 2017–2025 Rob Stoll
# Copyright (c) 2017–2025 Eric Pardyjak
# Copyright (c) 2017–2025 Pete Willemsen
#
# This file is part of UtahLSM.
#
# This software is free and is distributed under the MIT License.
# See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
#
"""Abstract base class for soil models in UtahLSM.

This module defines the data structures and interface for all soil physics
parameterizations. It provides:
    1.  A `SoilProperties` dataclass to hold soil parameter arrays.
    2.  The `Soil` abstract base class, which ensures that any concrete soil
        model implements the necessary methods.
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, fields
from typing import TypeVar, Union, overload

import numpy as np
from numpy.typing import NDArray

from .soil_type import SoilType
from ...util import constants as c
from ...util.io import logging_helper

ST = TypeVar('ST', bound='Soil')
logger = logging_helper.get_logger('SOIL')

@dataclass
class SoilProperties:
    """A container for soil property arrays.

    This dataclass holds the soil parameters for each layer in the soil
    column as NumPy arrays.

    Attributes:
        b: Clapp and Hornberger "b" exponent (unitless).
        psi_sat: Saturation moisture potential [m].
        porosity: Saturated soil moisture content (volumetric) [m^3/m^3].
        residual: Residual moisture content (volumetric) [m^3/m^3].
        K_sat: Saturated hydraulic conductivity [m/s].
        ci: Volumetric heat capacity [J/m^3-K].
    """
    b: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    psi_sat: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    porosity: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    residual: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    K_sat: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    ci: NDArray[np.float64] = field(default_factory=lambda: np.array([]))

class Soil(ABC):
    """Abstract base class for soil physics models.

    This class defines the standard interface for all soil models and acts
    as a factory. It also contains shared methods for calculating various
    soil-related quantities that are common across different parameterizations.

    Attributes:
        logger: A logger for this class.
        input: An `Input` object containing model configuration.
        properties: A `SoilProperties` object holding the soil parameters
            for each layer of the soil column.
    """
    def __init__(self, dataset_id: int, soil_type_array: NDArray[np.int_]):

        """Initializes the Soil model.

        This constructor populates the `properties` attribute by looking up
        the soil type for each layer from the input files and assembling the
        corresponding physical parameters into NumPy arrays.

        Args:
            dataset_id: An integer ID for the soil parameter dataset to use.
            soil_type_array: A NumPy array of soil type IDs for each layer.
        """
        self.logger: logging.Logger = logging_helper.get_logger('SOIL')

        dataset_names = {
            1: 'Clapp/Hornberger',
            2: 'Cosby et al',
            3: 'Rawls/Brakensiek'
        }
        self.logger.info('Using the %s dataset', dataset_names[dataset_id])

        # Create temporary lists to hold properties for each layer
        prop_lists = {f.name: [] for f in fields(SoilProperties)}

        # Loop to gather properties from the original SoilType objects
        for soil_type in soil_type_array:
            props = SoilType.get_properties(dataset_id, soil_type)
            for prop_name in prop_lists.keys():
                prop_lists[prop_name].append(getattr(props, prop_name))

        # Convert lists to arrays and store them in our dataclass
        self.properties = SoilProperties(
            **{name: np.array(values) for name, values in prop_lists.items()}
        )

    #--- Abstract Methods ---

    @overload
    def water_potential(self, soil_q: float, level: int = None) -> float: ...

    @overload
    def water_potential(
        self, soil_q: NDArray[np.float64], level: int = None
    ) -> NDArray[np.float64]: ...

    @abstractmethod
    def water_potential(
        self, soil_q: Union[float, NDArray[np.float64]], level: int = None
    ) -> Union[float, NDArray[np.float64]]:
        """Computes soil water potential using soil moisture content.

        This method must be implemented by subclasses to compute the water
        potential (matric potential) based on the soil moisture and soil
        water retention properties.

        Args:
            soil_q: Soil moisture content [m^3/m^3]. Can be a scalar value
                or an array of values.
            level: Optional specific soil layer index. If provided with a
                scalar soil_q, indicates which layer the moisture belongs to.

        Returns:
            Water potential [Pa]. Returns the same type as soil_q (float
            or NDArray).
        """
        raise NotImplementedError

    @overload
    def conductivity_moisture(
        self, soil_q: float, level: int = None
    ) -> float: ...

    @overload
    def conductivity_moisture(
        self, soil_q: NDArray[np.float64], level: int = None
    ) -> NDArray[np.float64]: ...

    @abstractmethod
    def conductivity_moisture(
        self, soil_q: Union[float, NDArray[np.float64]], level: int = None
    ) -> Union[float, NDArray[np.float64]]:
        """Computes soil moisture conductivity.

        This method must be implemented by subclasses to compute the water
        conductivity as a function of soil moisture using the soil's
        water retention and conductivity relationships.

        Args:
            soil_q: Soil moisture content [m^3/m^3]. Can be a scalar value
                or an array of values.
            level: Optional specific soil layer index. If provided with a
                scalar soil_q, indicates which layer the moisture belongs to.

        Returns:
            Moisture conductivity [m/s]. Returns the same type as soil_q
            (float or NDArray).
        """
        raise NotImplementedError

    @abstractmethod
    def diffusivity_moisture(
        self, soil_q: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Computes soil moisture diffusivity profile.

        This method must be implemented by subclasses to compute the moisture
        diffusivity as a function of soil moisture for all soil layers. The
        moisture diffusivity is used in the implicit diffusion solver for
        soil moisture transport.

        Args:
            soil_q: Soil moisture content profile [m^3/m^3]. Array with
                one element per soil layer.

        Returns:
            Moisture diffusivity profile [m^2/s]. Array with one element
            per soil layer.
        """
        raise NotImplementedError

    @abstractmethod
    def surface_water_content(self, psi_sfc: float) -> float:
        """Computes surface soil moisture from water potential.

        This method must be implemented by subclasses to compute the soil
        moisture at the surface given a water potential value. This is used
        to determine surface moisture content from the surface water
        potential solution.

        Args:
            psi_sfc: Water potential at the surface [Pa].

        Returns:
            Soil moisture content at the surface [m^3/m^3].
        """
        raise NotImplementedError

    # --- Validation Methods ---

    def _validate_moisture_bounds(
        self, soil_q: Union[float, NDArray[np.float64]], level: int = None
    ) -> None:
        """Validates that soil moisture is within physically possible bounds.

        Issues warnings for out-of-bounds values instead of raising errors,
        allowing simulation to continue with imperfect data for debugging.

        Args:
            soil_q: Soil moisture content [m^3/m^3].
            level: The specific soil layer index if soil_q is a scalar.
        """
        if level is not None:
            porosity = self.properties.porosity[level]
            if isinstance(soil_q, (int, float)):
                if soil_q < 0:
                    self.logger.warning(
                        'Layer %d: soil moisture %f is negative. '
                        'Moisture must be >= 0.', level, soil_q
                    )
                if soil_q > porosity:
                    self.logger.warning(
                        'Layer %d: soil moisture %f exceeds '
                        'porosity %f.', level, soil_q, porosity)
        else:
            porosity = self.properties.porosity
            if isinstance(soil_q, np.ndarray):
                invalid_neg = np.where(soil_q < 0)[0]
                if len(invalid_neg) > 0:
                    self.logger.warning(
                        'Layers %s have negative soil '
                        'moisture: %s. '
                        'Moisture must be >= 0.',
                        invalid_neg.tolist(), soil_q[invalid_neg].tolist())
                invalid_high = np.where(soil_q > porosity)[0]
                if len(invalid_high) > 0:
                    self.logger.warning(
                        'Layers %s have soil '
                        'moisture exceeding porosity: '
                        '%s > %s.',
                        invalid_high.tolist(), soil_q[invalid_high].tolist(),
                        porosity[invalid_high].tolist())

    # --- Shared Methods ---

    def heat_capacity(self, soil_q: NDArray[np.float64]) -> NDArray[np.float64]:
        """Computes the volumetric heat capacity of the soil.

        Args:
            soil_q: The soil moisture content for each layer [m^3/m^3].

        Returns:
            The volumetric heat capacity for each layer [J/m^3-K].
        """
        CI_W = c.water.SPECIFIC_HEAT
        CP_A = c.thermodynamic.SPECIFIC_HEAT

        porosity = self.properties.porosity
        Ci = self.properties.ci
        Ks = (1.-porosity)*Ci + soil_q*CI_W + (porosity-soil_q)*CP_A

        return Ks

    def surface_mixing_ratio(self, sfc_T: float, sfc_q: float,
                             atm_p: float) -> float:
        """Computes the specific humidity at the soil surface.

        Args:
            sfc_T: The surface temperature [K].
            sfc_q: The surface soil moisture content [m^3/m^3].
            atm_p: The atmospheric pressure [Pa].

        Returns:
            The specific humidity at the surface [kg/kg].
        """
        G  = c.physical.GRAVITY
        RV = c.thermodynamic.GAS_CONSTANT_VAPOR

        psi = self.water_potential(sfc_q, level=0)
        h = np.exp(G*psi/(RV*sfc_T))
        es = c.thermodynamic.ES_REF * np.exp(
            c.thermodynamic.TETENS_A * (sfc_T-c.air.TEMPERATURE_REF) /
            (sfc_T-c.thermodynamic.TETENS_B))
        hum_sat = c.thermodynamic.EPSILON*(es/(atm_p-0.378*es))
        hum_spec = h*hum_sat

        return hum_spec

    def conductivity_thermal(
            self, soil_q: NDArray[np.float64]) -> NDArray[np.float64]:
        """Computes the soil thermal conductivity.

        Args:
            soil_q: The soil moisture content for each layer [m^3/m^3].

        Returns:
            The thermal conductivity for each layer [W/m-K].
        """
        psi = self.water_potential(soil_q)
        pf = np.log10(np.abs(psi * 100) + 1e-9)
        conductivity = np.where(
            pf <= c.soil.CONDUCTIVITY_PF_THRESHOLD,
            c.soil.CONDUCTIVITY_COEFF * np.exp(-(pf + c.soil.CONDUCTIVITY_EXP)),
            c.soil.CONDUCTIVITY_MIN
        )
        return conductivity

    def diffusivity_thermal(
            self, soil_q: NDArray[np.float64]) -> NDArray[np.float64]:
        """Computes the soil thermal diffusivity.

        Args:
            soil_q: The soil moisture content for each layer [m^3/m^3].

        Returns:
            The thermal diffusivity for each layer [m^2/s].
        """
        heat_cap = self.heat_capacity(soil_q)
        conductivity = self.conductivity_thermal(soil_q)
        diffusivity = conductivity / heat_cap
        return diffusivity
