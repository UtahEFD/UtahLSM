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
"""Abstract base class for surface layer models in UtahLSM.

This module defines the interface for all surface layer physics schemes
used within the UtahLSM framework. It provides the `Surface` abstract base
class, which ensures that any concrete surface layer model implements the
necessary stability functions.
"""

import logging
from abc import ABC, abstractmethod
from typing import TypeVar

from ..._types import FloatOrArray
from ...util.io import logging_helper

ST = TypeVar('ST', bound='Surface')
logger = logging_helper.get_logger('SFC')

class Surface(ABC):
    """Abstract base class for surface layer models.

    This class defines the standard interface for surface layer calculations
    and acts as a factory for creating specific model instances. It cannot
    be instantiated directly.

    Attributes:
        logger: A logger for this class.
    """
    def __init__(self) -> None:
        """Initializes the Surface base class."""
        self.logger: logging.Logger = logging_helper.get_logger('SFC')

    #--- Abstract methods ---
    @abstractmethod
    def fm(self, z1: float, z0: float,
           obukhov_l: FloatOrArray
    ) -> FloatOrArray:
        """Computes the stability function for momentum.

        This is an abstract method that must be implemented by any concrete
        subclass.

        Args:
            z1: Upper height [m].
            z0: Lower height (roughness length) [m].
            obukhov_l: Obukhov length [m].

        Returns:
            The dimensionless stability correction factor for momentum.
        """
        raise NotImplementedError

    @abstractmethod
    def fh(self, z1: float, z0h: float,
           obukhov_l: FloatOrArray
    ) -> FloatOrArray:
        """Computes the stability function for heat.

        This is an abstract method that must be implemented by any concrete
        subclass.

        Args:
            z1: Upper height [m].
            z0h: Lower height (thermal roughness length) [m].
            obukhov_l: Obukhov length [m].

        Returns:
            The dimensionless stability correction factor for heat and scalars.
        """
        raise NotImplementedError
