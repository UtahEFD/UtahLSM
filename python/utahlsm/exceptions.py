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
"""Custom exception types for the UtahLSM model."""

class UtahLSMError(Exception):
    """Base class for all custom exceptions in the UtahLSM model."""
    pass

class NamelistError(UtahLSMError):
    """Raised for errors found in the namelist configuration."""
    pass