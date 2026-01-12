## Final Report: Covid in Canada SIRD

## Parameters Used
- **Population (POP)**: 38,940,000
- **Max Hospital Capacity (I_MAX)**: 62,000
- **Time Step (DT)**: 10 days

## For full data, see the `/data/populations.csv` file. (see sources below)

## Differentiation Results
- **Average dS/dt**: -4,029.32 people/day
- **Average dI/dt**: 79.14 people/day
- **Average dR/dt**: 3,901.70 people/day
- **Average dD/dt**: 48.48 people/day

## Interpolation Methods Used
- Lagrange Interpolation
- Linear Spline Interpolation
- Quadratic Spline Interpolation

## Integration Methods Used
- Trapezoidal Integration
- Simpson's Rule Integration

## Differential Equation Solvers Used
- Explicit Euler Method
- Implicit Euler Method

## Data Overview
- **Total Days Analyzed**: 750 days
- **Final Infected Count**: 59,360
- **Final Recovered Count**: 2,926,278
- **Final Deceased Count**: 36,359

- **Infection Rate (r)**: 2.1057e-09
- **Recovery Rate (a)**: 0.079802
- **Mortality Rate (b)**: 0.001326
- **Basic Reproduction Number (R0)**: 2.5955e-08

## Critical Indicators
- **tc (Hospital Saturation)**: 671.63 days
- **HIT (Collective Immunity Threshold)**: 412,082 people
- **Vaccine Doses Needed (70% HIT)**: 48,663,444

## Estimated Total Infections
- **Total Infections (Trapezoidal)**: 2,852,255
- **Total Infections (Simpson)**: 2,824,081

## Graphs
- See the `/graphs` directory for all generated plots.
