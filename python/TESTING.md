# UtahLSM Testing Guide

## Overview

This document provides guidance for using and extending the UtahLSM test suite. The test suite uses **pytest** as the testing framework and focuses on unit testing to validate individual components of the model.

## Test Organization

All tests are located in the `tests/` directory and organized by component:

```
tests/
├── conftest.py              # Shared fixtures and test utilities
├── __init__.py              # Test package initialization
├── test_solvers.py          # Numerical solvers (tridiagonal, root_brent)
├── test_soil_models.py      # Soil physics models
├── test_seb_solver.py       # Surface Energy Budget tests (future)
├── test_smb_solver.py       # Surface Moisture Budget tests (future)
└── test_data_models.py      # Configuration and state dataclasses
```

## Running Tests

### Run All Tests

```bash
# Run all tests with verbose output
pytest tests/ -v

# Run all tests and show coverage
pytest tests/ -v --cov=utahlsm
```

### Run Tests by Component

```bash
# Run only solver tests
pytest tests/test_solvers.py -v

# Run only data model tests
pytest tests/test_data_models.py -v

# Run only soil tests
pytest tests/test_soil_models.py -v
```

### Run Tests by Marker

```bash
# Run only solver tests (using pytest markers)
pytest -m solver

# Run only soil tests
pytest -m soil

# Run only data model tests
pytest -m datamodel
```

### Run Specific Tests

```bash
# Run a single test file
pytest tests/test_solvers.py::TestTridiagonal::test_single_element -v

# Run tests matching a pattern
pytest tests/ -k "monotonic" -v
```

## Test Fixtures

The `conftest.py` file provides reusable fixtures for testing. These fixtures initialize common test data structures:

### Atmospheric State Fixtures

- `atm_neutral`: Neutral atmospheric conditions (typical testing scenario)
- `atm_stable`: Stable atmospheric conditions (cold, calm)
- `atm_unstable`: Unstable atmospheric conditions (warm, weak wind)

### Soil State Fixtures

- `soil_state_simple`: Simple uniform soil state
- `soil_state_profile`: Realistic soil state with vertical gradients
- `soil_state_dry`: Dry soil state
- `soil_state_saturated`: Saturated soil state

### Surface State Fixtures

- `surface_state_base`: Basic surface state

### Configuration Fixtures

- `grid_config_5layers`: 5-layer soil grid
- `grid_config_10layers`: 10-layer soil grid (finer resolution)
- `surface_config`: Surface configuration with roughness parameters

### Helper Fixtures

- `create_tridiagonal_system`: Factory for generating test tridiagonal systems
- `create_root_function`: Factory for creating test functions with known roots

## Test Coverage by Component

### 1. Numerical Solvers (33 tests, 1 skipped)

**File:** `test_solvers.py`

Tests the core numerical solvers that are used throughout the model:

- **Tridiagonal Solver (Thomas Algorithm)**
  - Single and multiple element systems
  - Well-conditioned, diagonal-dominant, and ill-conditioned systems
  - Physical boundary conditions (heat diffusion analogy)
  - Error handling for zero on diagonal
  - Various matrix sizes (5-100 elements)
  - Symmetric matrices

- **Root-Finding Solver (Brent's Method)**
  - Linear, quadratic, cubic functions
  - Transcendental functions (sin, exp, tan)
  - Convergence with different tolerances
  - Maximum iteration limits
  - Bracketing errors
  - Energy balance analogy (physics-relevant test)
  - Oscillating and steep functions

**Key Physics Tests:**
- Heat diffusion with fixed boundary conditions
- Energy balance root-finding

### 2. Data Models (30 tests)

**File:** `test_data_models.py`

Tests the configuration and state dataclasses:

- **Atmospheric State**
  - Initialization with default and custom values
  - Reasonable temperature ranges
  - Physical constraints (wind speed, humidity, pressure)

- **Soil State**
  - Array consistency between temperature, moisture, and type
  - Realistic temperature and moisture ranges
  - Soil type tracking

- **Surface State**
  - Flux vectors and turbulence scales
  - Physically reasonable magnitude constraints

- **Grid Configuration**
  - Depth positivity
  - Monotonically increasing cumulative depths
  - Dimension consistency

- **Surface Configuration**
  - Roughness length relationships
  - Measurement height positivity
  - Albedo and emissivity bounds [0, 1]

- **State Constraints**
  - Temperature continuity between layers
  - Moisture bounds

### 3. Soil Models (18 tests, mostly skipped pending soil data loading)

**File:** `test_soil_models.py`

Tests soil physics models (Brooks-Corey, Campbell, Van Genuchten):

- **Water Potential**
  - Monotonic behavior with moisture
  - Reasonable values at saturation and residual
  - Physical bounds

- **Hydraulic Conductivity**
  - Matches K_sat at saturation
  - Decreases with dryness
  - Always positive

- **Thermal Properties**
  - Positive thermal conductivity and diffusivity
  - Increases with moisture

- **Parameter Constraints**
  - Residual < Porosity
  - Porosity ∈ [0, 1]
  - K_sat > 0
  - Heat capacity > 0

- **Numerical Stability**
  - No NaN near saturation
  - No Inf near residual moisture

## Key Testing Principles

### 1. Physics-First Approach

Tests focus on physical constraints and known relationships:

- Conservation laws (mass, energy)
- Monotonicity and sign constraints
- Boundary conditions
- Realistic value ranges

### 2. Analytical Solutions

Where possible, tests use analytical solutions:

- Tridiagonal solver: Tests against known solutions
- Root-finding: Tests with functions having known roots
- Energy balance: Tests with simplified but physics-relevant functions

### 3. Parametrized Tests

Many tests use `@pytest.mark.parametrize` to test multiple scenarios:

```python
@pytest.mark.parametrize("size", [5, 10, 20, 50])
def test_various_sizes(self, create_tridiagonal_system, size):
    a, b, c, r, x_expected = create_tridiagonal_system(size, 'well-conditioned')
    x = tridiagonal(a, b, c, r)
    assert_allclose(x, x_expected, rtol=1e-8)
```

### 4. Error Handling

Tests verify that the code properly handles invalid inputs:

- Zero on diagonal
- Unbracketed roots
- Physical constraint violations

## Understanding Test Results

### Test Output

```
tests/test_solvers.py::TestTridiagonal::test_single_element PASSED [ 50%]
```

- `tests/test_solvers.py`: Test file
- `TestTridiagonal`: Test class
- `test_single_element`: Test method
- `PASSED`: Result (PASSED, FAILED, SKIPPED)
- `[ 50%]`: Progress through test suite

### Common Issues

**SKIPPED Tests:**
- Tests marked as skipped typically require external data (e.g., soil properties)
- Check the skip message with `pytest -v` for details

**FAILED Tests:**
- Check the assertion error for details
- Use `pytest --tb=long` for full traceback
- Use `pytest -x` to stop on first failure for debugging

## Adding New Tests

### Test Structure

```python
@pytest.mark.soil
class TestNewFeature:
    """Test suite for new feature."""

    def test_specific_behavior(self, fixture_name):
        """Description of what is tested.

        Explanation of physics, constraints, or expected behavior.
        """
        # Setup
        input_data = fixture_name.some_method()

        # Execute
        result = calculate_something(input_data)

        # Assert
        assert result > 0.0, "Result should be positive"
        assert_allclose(result, expected, rtol=1e-8)
```

### Marker Usage

Use pytest markers to organize tests:

```python
@pytest.mark.solver      # For solver tests
@pytest.mark.soil        # For soil physics tests
@pytest.mark.seb         # For SEB solver tests
@pytest.mark.smb         # For SMB solver tests
@pytest.mark.datamodel   # For data model tests
@pytest.mark.unit        # For unit tests
@pytest.mark.integration # For integration tests
@pytest.mark.slow        # For slow-running tests
```

### Using Fixtures

```python
def test_with_fixtures(self, atm_neutral, soil_state_simple, grid_config_5layers):
    """Test that uses multiple fixtures."""
    # Fixtures are automatically instantiated and passed to test
    assert atm_neutral.wind_speed == 5.0
    assert soil_state_simple.temperature.shape == (5,)
```

## Integration with CI/CD

To add these tests to a CI/CD pipeline (GitHub Actions, etc.):

```yaml
# .github/workflows/tests.yml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
        with:
          python-version: '3.12'
      - run: pip install pytest numpy
      - run: pytest tests/ -v
```

## Best Practices

1. **One assertion per test method**: Each test should verify one behavior
2. **Clear test names**: Test names should describe what is tested
3. **Use fixtures**: Avoid code duplication with fixtures
4. **Test edge cases**: Test boundary conditions and extreme values
5. **Physical constraints**: Always verify that results are physically reasonable
6. **Document tests**: Include docstrings explaining the physics being tested

## Troubleshooting

### Import Errors

If you get import errors when running tests:

```bash
# Make sure the package is importable
python -c "import utahlsm; print(utahlsm.__file__)"

# Run from the python directory
cd /Users/jeremy.gibbs/Desktop/UtahLSM/python
pytest tests/ -v
```

### pytest Not Found

```bash
# Install pytest if needed
pip install pytest numpy pytest-cov
```

### Slow Tests

```bash
# Time each test
pytest tests/ -v --durations=10

# Skip slow tests
pytest tests/ -v -m "not slow"
```

## Future Testing Areas

The following components are ready for expansion:

1. **Surface Energy Budget Solver** (`test_seb_solver.py`)
   - Energy balance satisfaction
   - Flux magnitude constraints
   - Convergence of root-finding

2. **Surface Moisture Budget Solver** (`test_smb_solver.py`)
   - Water potential gradients
   - Moisture budget closure
   - Convergence of iterative solver

3. **Integration Tests**
   - Full simulation cycles with GABLS3 case
   - Output validation
   - Mass/energy conservation

4. **Regression Tests**
   - Reference output comparison
   - Performance benchmarks

## Questions or Issues?

For questions about the test suite:
1. Check docstrings in test files
2. Review the CLAUDE.md file for model architecture
3. Run `pytest --help` for command-line options
4. Check pytest documentation at https://docs.pytest.org/
