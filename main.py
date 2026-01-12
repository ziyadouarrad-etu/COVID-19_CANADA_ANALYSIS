import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.interpolations import lagrange_interpolation, linear_spline_interpolation, quadratic_spline_interpolation
from src.integrations import simpson_integration, trapezoidal_integration
from src.diff_equations import solve_sird_euler
from src.utils import bisection_method

# --- CONFIGURATION ---
POP = 38940000; I_MAX = 62000; DT = 10
os.makedirs('graphs', exist_ok=True)
df = pd.read_csv('data/populations.csv')
df['susceptible_can'] = POP - (df['infected_can'] + df['recovered_can'] + df['deceased_can'])

# --- 1. DIFFERENTIATION & PARAMETERS ---
diff = df.diff() / DT
a = (diff['recovered_can'] / df['infected_can']).mean()
b = (diff['deceased_can'] / df['infected_can']).mean()
r = -(diff['susceptible_can'] / (df['susceptible_can'] * df['infected_can'])).mean()
r0 = r / (a + b)

# --- 2. NON-LINEAR INDICATORS ---
# Hospital Saturation tc: Solve I(t) - I_MAX = 0
tc = bisection_method(df['day'], df['infected_can'], I_MAX, df['day'].min(), df['day'].max(), linear_spline_interpolation)
hit_val = (1 - (a + b) / (r * POP)) * POP

# --- 3. INTERPOLATION & VISUALIZATION ---
eval_x = np.linspace(0, 750, 1000)
for col, name in {'infected_can': 'Infected', 'deceased_can': 'Deceased', 'recovered_can': 'Recovered', 'susceptible_can': 'Susceptible'}.items():
    plt.figure(figsize=(10,6))
    plt.scatter(df['day'], df[col], color='red', label='Real Data')
    plt.plot(eval_x, linear_spline_interpolation(eval_x, df['day'], df[col]), 'b-', label='Linear Spline')
    plt.plot(eval_x, quadratic_spline_interpolation(eval_x, df['day'], df[col]), 'g--', label='Quadratic Spline')

    # Show Lagrange boundary instability (Runge)
    y_lag = [lagrange_interpolation(xi, df['day'], df[col]) for xi in eval_x]
    plt.plot(eval_x, y_lag, 'k:', alpha=0.3, label='Lagrange (Oscillations)')
    
    plt.ylim(df[col].min()-5000, df[col].max()+5000)
    plt.title(f'Analyse d\'Interpolation: {name}')
    plt.legend(); plt.grid(True); plt.savefig(f'graphs/{col}_analysis.png'); plt.close()

# --- 4. REGIONAL COMPARISON ---
df['dI_usa'] = df['infected_usa'].diff() / DT
plt.figure(figsize=(10,6))
plt.plot(df['day'], diff['infected_can'], label='Canada (vitesse)')
plt.plot(df['day'], df['dI_usa'], label='USA (vitesse)', linestyle='--')
plt.title('Comparaison de Vitesse de Propagation (dI/dt)')
plt.legend(); plt.savefig('graphs/regional_comparison.png'); plt.close()

# --- 5. Rate of change of Susceptible and Rate of change of Infected Plot in one plot ---
plt.figure(figsize=(10,6))
plt.plot(df['day'], diff['susceptible_can'], label='dS/dt (Susceptible Rate of Change)')
plt.plot(df['day'], diff['infected_can'], label='dI/dt (Infected Rate of Change)', linestyle='--')
plt.title('Rate of Change: Susceptible vs Infected')
plt.legend(); plt.savefig('graphs/rate_of_change.png'); plt.close()

# --- 6. INTEGRATION & TOTALS ---
# I_cumulative(t) = ∫ r*S(τ)*I(τ) dτ
df['new_inf_rate'] = r * df['susceptible_can'] * df['infected_can']
# Use Trapezoidal for Total Infections
total_infected_trap = trapezoidal_integration(
    df['new_inf_rate'], 
    initial_value=df['infected_can'].iloc[0], 
    start_idx=0, 
    end_idx=len(df)-1, 
    step=10
)
# Use Simpson for Total Infections
total_infected_simp = simpson_integration(
    df['new_inf_rate'],
    initial_value=df['infected_can'].iloc[0],
    start_idx=0,
    end_idx=len(df)-1,
    step=10
)

# R_cumulé(t) = ∫ a*I(τ) dτ
df['recov_rate'] = a * df['infected_can']
# Use Trapezoidal for Total Recovered
total_recovered_trap = trapezoidal_integration(
    df['recov_rate'], 
    initial_value=0, 
    start_idx=0, 
    end_idx=len(df)-1, 
    step=10
)
# Use Simpson for Total Recovered
total_recovered_simp = simpson_integration(
    df['recov_rate'],
    initial_value=0,
    start_idx=0,
    end_idx=len(df)-1,
    step=10
)

# --- 7. SIRD SOLUTION ---
# Solve SIRD using Explicit Euler
sird_results_explicit = solve_sird_euler(df, r, a, b, DT, method='explicit')
# Solve SIRD using Implicit Euler
sird_results_implicit = solve_sird_euler(df, r, a, b, DT, method='implicit')
# Plot SIRD results
plt.figure(figsize=(10,6))
plt.plot(sird_results_explicit['day'], sird_results_explicit['infected'], label='Infected (Explicit)', linestyle='-')
plt.plot(sird_results_implicit['day'], sird_results_implicit['infected'], label='Infected (Implicit)', linestyle='--')
plt.title('SIRD Model: Explicit vs Implicit Euler')
plt.legend(); plt.savefig('graphs/sird_comparison.png'); plt.close()

# --- 8. EXPORT RESULTS ---
with open("results.md", "w") as f:
    f.write(f"## Final Report: Covid in Canada SIRD\n\n")
    f.write(f"## Parameters Used\n- **Population (POP)**: {POP:,}\n- **Max Hospital Capacity (I_MAX)**: {I_MAX:,}\n- **Time Step (DT)**: {DT} days\n\n")
    f.write(f"## For full data, see the `/data/populations.csv` file. (see sources below)\n\n")
    f.write(f"## Differentiation Results\n- **Average dS/dt**: {diff['susceptible_can'].mean():,.2f} people/day\n- **Average dI/dt**: {diff['infected_can'].mean():,.2f} people/day\n- **Average dR/dt**: {diff['recovered_can'].mean():,.2f} people/day\n- **Average dD/dt**: {diff['deceased_can'].mean():,.2f} people/day\n\n")
    f.write(f"## Interpolation Methods Used\n- Lagrange Interpolation\n- Linear Spline Interpolation\n- Quadratic Spline Interpolation\n\n")
    f.write(f"## Integration Methods Used\n- Trapezoidal Integration\n- Simpson's Rule Integration\n\n")
    f.write(f"## Differential Equation Solvers Used\n- Explicit Euler Method\n- Implicit Euler Method\n\n")
    f.write(f"## Data Overview\n- **Total Days Analyzed**: {df['day'].iloc[-1]} days\n- **Final Infected Count**: {df['infected_can'].iloc[-1]:,.0f}\n- **Final Recovered Count**: {df['recovered_can'].iloc[-1]:,.0f}\n- **Final Deceased Count**: {df['deceased_can'].iloc[-1]:,.0f}\n\n")
    f.write(f"- **Infection Rate (r)**: {r:.4e}\n")
    f.write(f"- **Recovery Rate (a)**: {a:.6f}\n")
    f.write(f"- **Mortality Rate (b)**: {b:.6f}\n")
    f.write(f"- **Basic Reproduction Number (R0)**: {r0:.4e}\n\n")
    f.write(f"## Critical Indicators\n- **tc (Hospital Saturation)**: {tc:.2f} days\n")
    f.write(f"- **HIT (Collective Immunity Threshold)**: {hit_val:,.0f} people\n")
    f.write(f"- **Vaccine Doses Needed (70% HIT)**: {(POP*0.7 - df['recovered_can'].iloc[-1])*2:,.0f}\n")
    f.write(f"\n## Estimated Total Infections\n")
    f.write(f"- **Total Infections (Trapezoidal)**: {total_infected_trap:,.0f}\n")
    f.write(f"- **Total Infections (Simpson)**: {total_infected_simp:,.0f}\n")
    f.write(f"\n## Graphs\n")
    f.write(f"- See the `/graphs` directory for all generated plots.\n")

print("✅ All tasks completed. Results in results.md and /graphs.")