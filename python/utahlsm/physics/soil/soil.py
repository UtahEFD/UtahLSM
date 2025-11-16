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
    3.  A factory function (`get_model`) for creating instances of specific
        soil models.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, fields
import logging
import numpy as np
from numpy.typing import NDArray
from typing import TypeVar, Union, overload
from .soil_type import SoilType
from ...exceptions import NamelistError
from ...util import constants as c
from ...util.io import logging_helper

ST = TypeVar('ST', bound='Soil')
logger = logging_helper.get_logger("SOIL")

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
        self.logger: logging.Logger = logging_helper.get_logger("SOIL")
        
        dataset_names = {
            1: "Clapp/Hornberger",
            2: "Cosby et al",
            3: "Rawls/Brakensiek"
        }
        self.logger.info(f"Using the {dataset_names[dataset_id]} dataset")
        
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
    
    @staticmethod
    def get_model(key: int, dataset_id: int, 
                  soil_type_array: NDArray[np.int_])->ST:
        """Factory method to select and instantiate a soil model.
        
        Args:
            key: An integer ID for the soil model to use.
            dataset_id: An integer ID for the soil parameter dataset.
            soil_type_array: A NumPy array of soil type IDs for each layer.
        
        Returns:
            An instance of a concrete `Soil` subclass.
        
        Raises:
            KeyError: If the provided `key` is not a valid model ID.
        """
        # import soil sub-classes
        from .soil_brookscorey import BrooksCorey
        from .soil_campbell import Campbell
        from .soil_vangenuchten import VanGenuchten
        
        # dictionary to map keys to classes
        soil_models = {
            1: BrooksCorey,
            2: Campbell,
            3: VanGenuchten,
        }
        
        # return class or throw error
        try:
            # look up model class from dictionary
            return soil_models[key](dataset_id, soil_type_array)
        except KeyError as e:
            error_msg = f"{key} is an invalid soil model."
            logger.error("x"*62)
            logger.error(f"Namelist Error: {error_msg}")
            logger.error(f"Valid options are:")
            for k,v in soil_models.items():
                logger.error(f"\t{k} ({v.__name__})")
            logger.error("x"*62)
            raise NamelistError(error_msg)
    
    #--- Abstract Methods ---

    @overload
    def water_potential(self, soil_q: float, level: int = None) -> float: ...

    @overload
    def water_potential(self, soil_q: NDArray[np.float64], level: int = None) -> NDArray[np.float64]: ...

    @abstractmethod
    def water_potential(self, soil_q: Union[float, NDArray[np.float64]], level: int = None) -> Union[float, NDArray[np.float64]]:
        """Computes soil water potential. Must be implemented by subclasses."""
        raise NotImplementedError
    
    @overload
    def conductivity_moisture(self, soil_q: float, level: int = None) -> float: ...

    @overload
    def conductivity_moisture(self, soil_q: NDArray[np.float64], level: int = None) -> NDArray[np.float64]: ...

    @abstractmethod
    def conductivity_moisture(self, soil_q: Union[float, NDArray[np.float64]], level: int = None) -> Union[float, NDArray[np.float64]]:
        """Computes soil moisture conductivity. Must be implemented by subclasses."""
        raise NotImplementedError
    
    @abstractmethod
    def diffusivity_moisture(self, soil_q: NDArray[np.float64]) -> NDArray[np.float64]:
        """Computes soil moisture diffusivity. Must be implemented by subclasses."""
        raise NotImplementedError
    
    @abstractmethod
    def surface_water_content(self, psi_sfc: float) -> float:
        """Computes sfc water content from potential. Must be implemented by subclasses."""
        raise NotImplementedError
    
    # --- Validation Methods ---

    def _validate_moisture_bounds(self, soil_q: Union[float, NDArray[np.float64]], level: int = None) -> None:
        """Validates that soil moisture is within physically possible bounds.

        Issues warnings for out-of-bounds values instead of raising errors,
        allowing simulation to continue with imperfect data for debugging.

        Args:
            soil_q: Soil moisture content [m^3/m^3].
            level: The specific soil layer index if soil_q is a scalar.
        """
        if level is not None:
            residual = self.properties.residual[level]
            porosity = self.properties.porosity[level]
            if isinstance(soil_q, (int, float)):
                if soil_q < 0:
                    self.logger.warning(
                        f"Layer {level}: soil moisture {soil_q} is negative. "
                        f"Moisture must be >= 0."
                    )
                if soil_q > porosity:
                    self.logger.warning(
                        f"Layer {level}: soil moisture {soil_q} exceeds porosity {porosity}."
                    )
        else:
            residual = self.properties.residual
            porosity = self.properties.porosity
            if isinstance(soil_q, np.ndarray):
                invalid_neg = np.where(soil_q < 0)[0]
                if len(invalid_neg) > 0:
                    self.logger.warning(
                        f"Layers {invalid_neg.tolist()} have negative soil moisture: "
                        f"{soil_q[invalid_neg].tolist()}. Moisture must be >= 0."
                    )
                invalid_high = np.where(soil_q > porosity)[0]
                if len(invalid_high) > 0:
                    self.logger.warning(
                        f"Layers {invalid_high.tolist()} have soil moisture exceeding porosity: "
                        f"{soil_q[invalid_high].tolist()} > {porosity[invalid_high].tolist()}."
                    )

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
        
    def surface_mixing_ratio(self, sfc_T: float, sfc_q: float, atm_p: float) -> float:
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
        es = c.thermodynamic.ES_REF*np.exp(c.thermodynamic.TETENS_A*(sfc_T-c.air.TEMPERATURE_REF)/(sfc_T-c.thermodynamic.TETENS_B))
        hum_sat = c.thermodynamic.EPSILON*(es/(atm_p-0.378*es))
        hum_spec = h*hum_sat
         
        return hum_spec
    
    def conductivity_thermal(self, soil_q: NDArray[np.float64]) -> NDArray[np.float64]:
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
    
    def diffusivity_thermal(self, soil_q: NDArray[np.float64]) -> NDArray[np.float64]:
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