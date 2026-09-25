import numpy as np
import spiceypy as spice


# Load kernels
spice.furnsh("kernels.tm")

def get_target_coverage(target):
    code = spice.bodn2c(target)

    for i in range(spice.ktotal("SPK")):
        file, _, _, _= spice.kdata(i, "SPK")

        ids = spice.spkobj(file)

        if code in ids:
            cover = spice.spkcov(file, code)
            return cover[0], cover[-1]

    return None


# Supported planets mapping to SPICE target and IAU frame
VALID_PLANETS = {
    "MERCURY": ("MERCURY", "IAU_MERCURY"),
    "VENUS":   ("VENUS", "IAU_VENUS"),
    "EARTH":   ("EARTH", "IAU_EARTH"),
    "MOON":    ("MOON", "IAU_MOON"),
    "MARS":    ("MARS BARYCENTER", "IAU_MARS"),
    "JUPITER": ("JUPITER BARYCENTER", "IAU_JUPITER"),
    "SATURN":  ("SATURN BARYCENTER", "IAU_SATURN"),
    "URANUS":  ("URANUS BARYCENTER", "IAU_URANUS"),
    "NEPTUNE": ("NEPTUNE BARYCENTER", "IAU_NEPTUNE"),
    "PLUTO":   ("PLUTO BARYCENTER", "IAU_PLUTO"),
}

# 1. Asks user for planet and loops until a valid planet is given
while True:
    planet = input("Enter planet (Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto): ").strip() or "Earth"
    planet_upper = planet.upper()
    if planet_upper in VALID_PLANETS:
        target, frame = VALID_PLANETS[planet_upper]
        break
    print(f"'{planet}' is not recognized. Please choose from: {', '.join(p.capitalize() for p in VALID_PLANETS.keys())}\n")


coverage = get_target_coverage(target)

if coverage:
    start_et, end_et = coverage
    print("\nValid dates for this object:")
    print("  Earliest:", spice.et2utc(start_et, "C", 3))
    print("  Latest:  ", spice.et2utc(end_et, "C", 3))


# 2. Asks user for date and loops until a valid calendar date is given
# Ephemeris time is a clock that is strictly uninterrupted unlike UTC or MST
while True:
    date_str = input("Enter date (2026-09-09 12:00:00 UTC): ").strip() or "2026-09-09 12:00:00 UTC"
    try:
        et = spice.str2et(date_str)
        break
    except spice.utils.exceptions.SpiceyError:
        print(f"Invalid date '{date_str}'. Please enter a valid calendar date (e.g. 2026-09-09 12:00:00 UTC).\n")

# -----------------------------------------------------------------------------
# 3. Position & Velocity Calculation
# -----------------------------------------------------------------------------
# spice.spkezr returns the 6-element state vector [x, y, z, vx, vy, vz]
# relative to the Sun in the non-rotating J2000 inertial frame.
state, lt = spice.spkezr(target, et, "J2000", "NONE", "SUN")
position = state[:3]  # [x, y, z] Cartesian position in kilometers
dist_km = np.linalg.norm(position)  # Euclidean distance: sqrt(x^2 + y^2 + z^2)

print(f"\n--- {planet.capitalize()} on {date_str} ---")
# Physical position: Coordinates from the center of the Sun (0, 0, 0)
print(f"Position (x, y, z) [km]: {position}")
# Distance from Sun in km and Astronomical Units (1 AU ≈ 149.6 million km)
print(f"Distance to Sun:        {dist_km:,.1f} km ({dist_km / 1.496e8:.3f} AU)")

# -----------------------------------------------------------------------------
# 4. Planetary Orientation & Rotation
# -----------------------------------------------------------------------------
# pxform calculates the 3x3 rotation matrix to transform vectors from
# space (J2000) into coordinates on the planet's rotating surface (IAU frame).
rot_matrix = spice.pxform("J2000", frame, et)

# The planet's spin axis (North Pole) is [0, 0, 1] in its body-fixed frame.
# Transforming this vector into J2000 space gives the direction of its pole:
pole = rot_matrix.T @ np.array([0.0, 0.0, 1.0])
_, ra, dec = spice.recrad(pole)  # Converts (x, y, z) into Right Ascension & Declination

print(f"\nOrientation Frame:      {frame}")
# North Pole in J2000: Shows which point in the night sky the planet's spin axis points to.
# Dec (Declination): Angle above Earth's celestial equator (-90° to +90°).
# RA (Right Ascension): Celestial longitude along the celestial equator (0° to 360°).
print(f"North Pole in J2000:    RA = {np.degrees(ra):.2f}°, Dec = {np.degrees(dec):.2f}°")

# The rows of the rotation matrix represent the planet's axes in space at this second:
#   X = Direction of Prime Meridian (0° Lat, 0° Lon) in space (spins 360° each day)
#   Y = Direction of 90° East Longitude on the equator
#   Z = Direction of the North Pole spin axis
print(f"Rotation Matrix (J2000 -> {frame}):")
print(f"  X = [{rot_matrix[0, 0]:.6f}, {rot_matrix[0, 1]:.6f}, {rot_matrix[0, 2]:.6f}]")
print(f"  Y = [{rot_matrix[1, 0]:.6f}, {rot_matrix[1, 1]:.6f}, {rot_matrix[1, 2]:.6f}]")
print(f"  Z = [{rot_matrix[2, 0]:.6f}, {rot_matrix[2, 1]:.6f}, {rot_matrix[2, 2]:.6f}]")

# Unload all kernels from memory
spice.kclear()