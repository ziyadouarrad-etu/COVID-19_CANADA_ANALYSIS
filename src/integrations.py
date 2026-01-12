import numpy as np

def trapezoidal_integration(df_column, initial_value, start_idx, end_idx, step=10):
    # Formula: I_total = ∫(dI/dt) dt ≈ Σ [ (f(xi) + f(xi+1)) / 2 ] * h
    total = initial_value
    # df_column should be the derivative dI/dt
    for i in range(start_idx, end_idx):
        deriv = df_column.iloc[i]
        next_deriv = df_column.iloc[i + 1]
        # Calculation happening inside: Area of the trapezoid
        total += ((deriv + next_deriv) / 2) * step
    return total

def simpson_integration(df_column, initial_value, start_idx, end_idx, step=10):
    # Formula: (h/3) * [f(x0) + 4f(x1) + f(x2)]
    # Ensure we have an even number of intervals
    if (end_idx - start_idx) % 2 != 0:
        end_idx -= 1
    
    total = initial_value
    for i in range(start_idx, end_idx, 2):
        d0 = df_column.iloc[i]
        d1 = df_column.iloc[i + 1]
        d2 = df_column.iloc[i + 2]
        # Calculation happening inside: Simpson's weighted sum
        total += (step / 3) * (d0 + 4 * d1 + d2)
    return total