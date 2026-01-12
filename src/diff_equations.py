import numpy as np
import pandas as pd

def solve_sird_euler(df, r, a, b, dt=10, method='explicit'):
    # Solves the SIRD system using Explicit or Implicit Euler methods.

    n = len(df)
    # Initialize results containers
    S = np.zeros(n); I = np.zeros(n); R = np.zeros(n); D = np.zeros(n)
    
    # Set initial conditions from provided data
    S[0] = df['susceptible_can'].iloc[0]
    I[0] = df['infected_can'].iloc[0]
    R[0] = df['recovered_can'].iloc[0]
    D[0] = df['deceased_can'].iloc[0]

    for t in range(1, n):
        if method == 'explicit':
            # Explicit Euler Formulas
            S[t] = S[t-1] - dt * r * S[t-1] * I[t-1]
            I[t] = I[t-1] + dt * (r * S[t-1] * I[t-1] - (a + b) * I[t-1])
            R[t] = R[t-1] + dt * a * I[t-1]
            D[t] = D[t-1] + dt * b * I[t-1]
        else:
            # Implicit Euler Formulas
            S[t] = S[t-1] / (1 + dt * r * I[t-1])
            I[t] = I[t-1] / (1 - dt * r * S[t] + dt * (a + b))
            R[t] = R[t-1] + dt * a * I[t]
            D[t] = D[t-1] + dt * b * I[t]
    # Return results as DataFrame containing day, susceptible, infected, recovered, deceased
    results = pd.DataFrame({
        'day': df['day'],
        'susceptible': S,
        'infected': I,
        'recovered': R,
        'deceased': D
    })
    return results