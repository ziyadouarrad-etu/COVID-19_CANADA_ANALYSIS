def bisection_method(x_pts, y_pts, target, low, high, function, tol=1):
    """
    Implements the Dichotomy method to find when a population hits a 'target'.
    Used for finding critical intervention times (tc).
    function: interpolation function to evaluate (e.g., linear_spline_interpolation)
    """
    while abs(high - low) > tol:
        mid = (low + high) / 2
        val = function([mid], x_pts, y_pts)[0] - target
        
        # Check sign change
        if (function([low], x_pts, y_pts)[0] - target) * val < 0:
            high = mid
        else:
            low = mid
    return mid