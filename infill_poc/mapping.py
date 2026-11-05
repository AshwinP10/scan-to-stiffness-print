import numpy as np

K_MIN, K_MAX = 0.5, 3.0  # N/mm, placeholder probe range

# Capped below 100% because PrusaSlicer turns 100% sparse infill into solid infill,
# which is a different extrusion type and would drop out of the verification.
INFILL_MIN, INFILL_MAX = 10.0, 90.0


def stiffness_to_infill(k, k_min=K_MIN, k_max=K_MAX, infill_min=INFILL_MIN, infill_max=INFILL_MAX):
    """Placeholder linear stiffness (N/mm) -> infill (%) map; swap for the measured calibration curve."""
    t = np.clip((np.asarray(k, dtype=float) - k_min) / (k_max - k_min), 0.0, 1.0)
    return infill_min + t * (infill_max - infill_min)
