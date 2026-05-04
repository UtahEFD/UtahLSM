#
# UtahLSM
#
# Copyright (c) 2017–2026 Jeremy A. Gibbs
# Copyright (c) 2017–2026 Rob Stoll
# Copyright (c) 2017–2026 Eric Pardyjak
# Copyright (c) 2017–2026 Pete Willemsen
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
from typing import Any, TypeVar

import numpy as np
from numpy.typing import NDArray

from ..._types import FloatOrArray
from ...exceptions import NamelistError, SolverError
from ...util import constants as c
from ...util.io import logging_helper
from .. import thermo

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
    def __init__(
        self,
        properties_dict: dict[str, dict[str, Any]],
        soil_type_names: list[str],
        dataset_name: str = 'custom'
    ) -> None:
        """Initializes the Soil model.

        This constructor populates the `properties` attribute from a pre-loaded
        properties dictionary. The dictionary maps soil type names to their
        property dictionaries.

        Args:
            properties_dict: Dictionary mapping soil type names (lowercase strings)
                to property dicts. Each property dict must contain keys:
                'b', 'psi_sat', 'porosity', 'residual', 'K_sat', 'ci'.
            soil_type_names: List of soil type names for each layer (e.g.,
                ['sand', 'loam', 'clay']).
            dataset_name: Human-readable name of the dataset being used
                (for logging). Defaults to 'custom'.

        Raises:
            NamelistError: If a soil type is not found in properties_dict.
        """
        self.logger: logging.Logger = logging_helper.get_logger('SOIL')

        self.logger.info('Using soil property dataset: %s', dataset_name)

        # Create temporary lists to hold properties for each layer
        prop_lists: dict[str, list[float]] = {
            f.name: [] for f in fields(SoilProperties)
        }

        # Loop to gather properties from the pre-loaded dictionary
        for soil_type_name in soil_type_names:
            soil_type_lower = soil_type_name.lower()

            if soil_type_lower not in properties_dict:
                available = ', '.join(sorted(properties_dict.keys()))
                raise NamelistError(
                    f"Soil type '{soil_type_name}' not found in property set "
                    f"'{dataset_name}'. Available types: {available}"
                )

            props = properties_dict[soil_type_lower]
            for prop_name in prop_lists.keys():
                if prop_name not in props:
                    raise NamelistError(
                        f"Missing property '{prop_name}' for soil type "
                        f"'{soil_type_name}' in dataset '{dataset_name}'"
                    )
                prop_lists[prop_name].append(props[prop_name])

        # Convert lists to arrays and store them in our dataclass
        self.properties = SoilProperties(
            **{name: np.array(values) for name, values in prop_lists.items()}
        )

        # Derived per-layer quantities from the retention curve.
        # Field capacity and wilting point correspond to standard matric
        # potentials (-3.3 m ≈ -33 kPa and -150 m ≈ -1500 kPa) inverted
        # through the soil model's own ψ(θ) relation — they therefore stay
        # consistent with the dataset's retention parameters and with
        # whichever soil model is active.
        self.theta_fc: NDArray[np.float64] = self._theta_at_potential(
            c.soil.PSI_FIELD_CAPACITY
        )
        self.theta_wilt: NDArray[np.float64] = self._theta_at_potential(
            c.soil.PSI_WILTING_POINT
        )

    def _theta_at_potential(
        self, psi_target: float
    ) -> NDArray[np.float64]:
        """Returns soil moisture at a specified matric potential per layer.

        Inverts ψ(θ) for each soil layer via bisection on the interval
        [residual, porosity]. Called once at init to populate θ_fc and
        θ_wilt; not performance-critical.

        Args:
            psi_target: Target matric potential [m] (negative).

        Returns:
            Per-layer soil moisture [m^3/m^3] of shape (nz,) that yields
            ψ(θ) ≈ psi_target in each layer's retention curve.
        """
        nz = self.properties.b.size
        out = np.empty(nz)
        for i in range(nz):
            lo = self.properties.residual[i] + 1e-6
            hi = self.properties.porosity[i] - 1e-6
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                psi_mid = float(self.water_potential(mid, level=i))
                # ψ is monotonically increasing in θ (less negative as θ↑)
                if psi_mid < psi_target:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1e-8:
                    break
            out[i] = 0.5 * (lo + hi)
        return out

    #--- Abstract Methods ---

    @abstractmethod
    def water_potential(
        self, soil_q: FloatOrArray, level: int | None = None
    ) -> FloatOrArray:
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

    @abstractmethod
    def water_content(
        self, psi: FloatOrArray, level: int | None = None
    ) -> FloatOrArray:
        """Computes soil moisture content from water potential.

        This is the inverse of ``water_potential`` and is required by the
        mixed-form Richards solver, which iterates in pressure head while
        storing moisture as the prognostic state.

        Args:
            psi: Soil water potential [m]. Can be a scalar value or an array
                of values.
            level: Optional specific soil layer index. If provided with a
                scalar ``psi``, indicates which layer the potential belongs to.

        Returns:
            Soil moisture content [m^3/m^3]. Returns the same type as ``psi``.
        """
        raise NotImplementedError

    @abstractmethod
    def moisture_capacity(
        self, psi: FloatOrArray, level: int | None = None
    ) -> FloatOrArray:
        """Computes the specific moisture capacity dθ/dψ.

        Args:
            psi: Soil water potential [m]. Can be a scalar value or an array
                of values.
            level: Optional specific soil layer index. If provided with a
                scalar ``psi``, indicates which layer the potential belongs to.

        Returns:
            Specific moisture capacity [m^3/m^3 per m head]. Returns the same
            type as ``psi``.
        """
        raise NotImplementedError

    @abstractmethod
    def conductivity_moisture(
        self, soil_q: FloatOrArray, level: int | None = None
    ) -> FloatOrArray:
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
    def conductivity_gradient(
        self, soil_q: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Computes the linearized hydraulic conductivity gradient dK/dθ.

        Returns the secant linearization K(Se)/Se/(φ-θ_r) appropriate for
        each soil model, used in the gravity drainage term of the moisture
        diffusion solver.

        Args:
            soil_q: Soil moisture content profile [m^3/m^3]. Array with
                one element per soil layer.

        Returns:
            Linearized dK/dθ profile [m/s]. Array with one element per
            soil layer.
        """
        raise NotImplementedError

    def surface_water_content(self, psi_sfc: float) -> float:
        """Computes surface soil moisture from water potential.

        This is a convenience wrapper for the top soil layer's inverse
        retention curve.

        Args:
            psi_sfc: Water potential at the surface [Pa].

        Returns:
            Soil moisture content at the surface [m^3/m^3].
        """
        return float(self.water_content(psi_sfc, level=0))

    # --- Validation Methods ---

    def enforce_moisture_bounds(
        self, moisture: NDArray[np.float64], tol: float = 1e-8
    ) -> None:
        """Enforces physical moisture bounds after a solve step.

        Raises SolverError if any value lies outside [residual - tol,
        porosity + tol]. Clips and warns for soft overshoots within tol.

        Args:
            moisture: Soil moisture profile [m^3/m^3], modified in-place.
            tol: Tolerance for soft clipping before a hard error is raised.
        """
        moisture_arr = np.asarray(moisture, dtype=float)
        _2d = moisture_arr.ndim == 2
        residual = self.properties.residual[:, None] if _2d else self.properties.residual
        porosity = self.properties.porosity[:, None] if _2d else self.properties.porosity
        tol = max(tol, 1e-8)

        below_hard = moisture_arr < (residual - tol)
        above_hard = moisture_arr > (porosity + tol)
        hard_mask = below_hard | above_hard
        if np.any(hard_mask):
            min_delta = float(np.min(moisture_arr - residual))
            max_delta = float(np.max(moisture_arr - porosity))
            bad_idx = np.argwhere(hard_mask)
            examples = ", ".join(
                f"{tuple(int(i) for i in idx)}"
                f"={float(moisture_arr[tuple(idx)]):.6f}"
                for idx in bad_idx[:5]
            )
            raise SolverError(
                'Soil moisture left physical bounds after the mixed moisture '
                f'solve: {int(np.sum(hard_mask))} cells outside '
                f'[residual, porosity] by more than tol={tol:.1e}. '
                f'Min(theta-residual)={min_delta:.3e}, '
                f'Max(theta-porosity)={max_delta:.3e}. '
                f'Examples: {examples}'
            )

        clip_mask = (moisture_arr < residual) | (moisture_arr > porosity)
        if np.any(clip_mask):
            self.logger.warning(
                'Clipping %d soil moisture values to [residual, porosity] '
                'after the mixed moisture solve (tol=%.1e).',
                int(np.sum(clip_mask)),
                tol,
            )
            np.clip(moisture, residual, porosity, out=moisture)

    # --- Shared Methods ---

    def heat_capacity(self, soil_q: NDArray[np.float64]) -> NDArray[np.float64]:
        """Computes the volumetric heat capacity of the soil.

        Args:
            soil_q: The soil moisture content for each layer [m^3/m^3].

        Returns:
            The volumetric heat capacity for each layer [J/m^3-K].
        """
        CI_W = c.water.VOLUMETRIC_HEAT_CAPACITY
        CI_A = c.air.DENSITY_REF * c.thermodynamic.SPECIFIC_HEAT

        _2d = soil_q.ndim == 2
        porosity = self.properties.porosity[:, None] if _2d else self.properties.porosity
        Ci = self.properties.ci[:, None] if _2d else self.properties.ci
        Ks = (1.-porosity)*Ci + soil_q*CI_W + (porosity-soil_q)*CI_A

        return Ks

    def surface_specific_humidity(
        self,
        sfc_T: FloatOrArray,
        sfc_theta: FloatOrArray,
        atm_p: FloatOrArray
    ) -> FloatOrArray:
        """Computes the specific humidity at the soil surface.

        Args:
            sfc_T: The surface temperature [K].
            sfc_theta: The surface volumetric soil moisture [m^3/m^3].
            atm_p: The atmospheric pressure [Pa].

        Returns:
            The specific humidity at the surface [kg/kg].
        """
        G  = c.physical.GRAVITY
        RV = c.thermodynamic.GAS_CONSTANT_VAPOR

        psi = self.water_potential(sfc_theta, level=0)
        h = np.exp(G*psi/(RV*sfc_T))
        q_sat = thermo.saturation_specific_humidity(sfc_T, atm_p)

        return h * q_sat

    def conductivity_thermal(
            self, soil_q: NDArray[np.float64]) -> NDArray[np.float64]:
        """Computes the soil thermal conductivity.

        Args:
            soil_q: The soil moisture content for each layer [m^3/m^3].

        Returns:
            The thermal conductivity for each layer [W/m-K].
        """
        psi = np.asarray(self.water_potential(soil_q))
        pf = np.log10(np.abs(psi * 100) + 1e-9)
        conductivity: NDArray[np.float64] = np.where(
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
