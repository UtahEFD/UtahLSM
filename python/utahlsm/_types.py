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
"""Shared type aliases for UtahLSM."""

from collections.abc import Sequence
from typing import Union

import numpy as np
from numpy.typing import NDArray

FloatOrArray = Union[float, NDArray[np.float64]]
FloatOrArrayLike = Union[float, Sequence[float], NDArray[np.float64]]
