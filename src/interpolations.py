import pandas as pd
import numpy as np

def lagrange_interpolation(x_eval, x_pts, y_pts):
    n = len(x_pts)
    result = 0.0
    for i in range(n):
        xi, yi = x_pts[i], y_pts[i]
        term = yi
        for j in range(n):
            if j != i:
                xj = x_pts[j]
                term *= (x_eval - xj) / (xi - xj)
        result += term
    return result

def linear_spline_interpolation(x_eval, x_pts, y_pts):
    # Connects data points (xi, yi) by line segments.
    # Calculates A (slope) and B (intercept) for each segment.
    n = len(x_pts)
    
    # Initialize slope (A) and intercept (B) arrays
    A = np.zeros(n-1)
    B = np.zeros(n-1)
    
    # Calculate spline coefficients for each interval
    for i in range(n-1):
        # A[i] = (y_next - y) / (x_next - x)
        A[i] = (y_pts[i+1] - y_pts[i]) / (x_pts[i+1] - x_pts[i])
        # B[i] = y - A * x
        B[i] = y_pts[i] - A[i] * x_pts[i]
    
    results = []
    for x in x_eval:
        # Find the correct segment for the current x value
        # searchsorted efficiently finds the interval index
        idx = np.searchsorted(x_pts, x) - 1
        idx = np.clip(idx, 0, n-2)
        
        # Calculate y = Ax + B
        y = A[idx] * x + B[idx]
        results.append(y)
        
    return results


def quadratic_spline_interpolation(x_eval, x_pts, y_pts):
    # Interpolates points using degree-2 polynomials with continuity conditions.
    n = len(x_pts)
    
    A = np.zeros(n-1)
    B = np.zeros(n-1)
    C = np.zeros(n-1)
    
    # 1. Initialize first spline (Assuming second derivative is zero at start)
    A[0] = 0
    B[0] = (y_pts[1] - y_pts[0]) / (x_pts[1] - x_pts[0])
    C[0] = y_pts[0]
    
    # 2. Calculate subsequent splines recursively
    for i in range(1, n-1):
        dx = x_pts[i+1] - x_pts[i]
        dy = y_pts[i+1] - y_pts[i]
        
        # Continuity of first derivative and value at node i
        # A[i] = ((y_next - y) / dx^2) - ((2*A_prev*x + B_prev) / dx)
        A[i] = (dy / (dx**2)) - (2*A[i-1]*x_pts[i] + B[i-1]) / dx
        B[i] = 2*A[i-1]*x_pts[i] + B[i-1] - 2*A[i]*x_pts[i]
        C[i] = y_pts[i] - A[i]*(x_pts[i]**2) - B[i]*x_pts[i]
    
    # 3. Vectorized Evaluation
    results = []
    for x in x_eval:
        # Efficiently find segment index
        idx = np.searchsorted(x_pts, x) - 1
        idx = np.clip(idx, 0, n-2)
        
        # Quadratic formula: y = Ax^2 + Bx + C
        y = A[idx]*(x**2) + B[idx]*x + C[idx]
        results.append(y)
        
    return results
