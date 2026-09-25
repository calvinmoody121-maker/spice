import numpy as np
import spiceypy as spice


# =============================================================================
# STEP 1: Load SPICE Kernels into Memory
# =============================================================================
# spice.furnsh() loads kernel files into the internal SPICE "kernel pool".
# "kernels.tm" is a meta-kernel that lists all individual kernels to load:
#   - naif0012.tls: Leapseconds Kernel (LSK) for UTC <-> Ephemeris Time conversion
#   - de440.bsp:    Binary Ephemeris (SPK) containing planet position & velocity data
#   - pck00011.tpc: Text PCK containing planet sizes and rotation/orientation models
#   - *.bpc / *.tf: High-precision Earth and Moon orientation kernels
spice.furnsh("kernels.tm")


# =============================================================================
# Helper: Check the Valid Date Range Stored Inside the SPK Ephemeris Kernel
# =============================================================================
def get_target_coverage(target):
    # 1. spice.bodn2c: "Body Name to Code" translates a string like "EARTH" 
    #    into its NAIF integer ID (e.g., Earth = 399, Sun = 10, Mars Barycenter = 4)
    code = spice.bodn2c(target)

    # 2. spice.ktotal: Counts how many kernels of type "SPK" (orbit data) are loaded
    for i in range(spice.ktotal("SPK")):
        # 3. spice.kdata: Retrieves the filename of the loaded SPK kernel (e.g. "de440.bsp")
        file, _, _, _ = spice.kdata(i, "SPK")

        # 4. spice.spkobj: Scans the SPK file and returns a list of all body IDs it contains
        ids = spice.spkobj(file)

        # 5. If the target planet exists in this file, query its coverage interval
        if code in ids:
            # spice.spkcov: Returns a SPICE window/array containing start & end times (in ET)
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

# =============================================================================
# STEP 2: Ask User for Target Planet & Inspect its Coverage
# =============================================================================
while True:
    planet = input("Enter planet (Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto): ").strip() or "Earth"
    planet_upper = planet.upper()
    if planet_upper in VALID_PLANETS:
        target, frame = VALID_PLANETS[planet_upper]
        break
    print(f"'{planet}' is not recognized. Please choose from: {', '.join(p.capitalize() for p in VALID_PLANETS.keys())}\n")

# Query the kernel file to find the earliest and latest available dates
coverage = get_target_coverage(target)
if coverage:
    start_et, end_et = coverage
    print("\nValid dates for this object:")
    # spice.et2utc converts internal Ephemeris Time seconds back into a calendar string:
    #   "C" = Calendar format, 3 = decimal places for seconds
    print("  Earliest:", spice.et2utc(start_et, "C", 3))
    print("  Latest:  ", spice.et2utc(end_et, "C", 3))


# =============================================================================
# STEP 3: Ask User for Date & Convert to Ephemeris Time (ET)
# =============================================================================
# Ephemeris Time (ET/TDB) is a continuous, uniform clock used by orbital mechanics,
# measuring seconds past J2000 (Jan 1, 2000 12:00:00 TDB).
while True:
    date_str = input("Enter date (2026-09-09 12:00:00 UTC): ").strip() or "2026-09-09 12:00:00 UTC"
    try:
        # spice.str2et uses the Leapseconds Kernel (naif0012.tls) to convert
        # the human calendar string into Ephemeris Time (ET seconds past J2000)
        et = spice.str2et(date_str)
        break
    except spice.utils.exceptions.SpiceyError:
        print(f"Invalid date '{date_str}'. Please enter a valid calendar date (e.g. 2026-09-09 12:00:00 UTC).\n")

# =============================================================================
# STEP 4: Query Ephemeris Kernel for 3D Position & Velocity
# =============================================================================
# spice.spkezr queries the binary ephemeris kernel (de440.bsp) and returns
# the 6-element state vector [x, y, z, vx, vy, vz] relative to the Sun.
state, lt = spice.spkezr(target, et, "J2000", "NONE", "SUN")
position = state[:3]  # [x, y, z] Cartesian position in kilometers
dist_km = np.linalg.norm(position)  # Euclidean distance: sqrt(x^2 + y^2 + z^2)

print(f"\n--- {planet.capitalize()} on {date_str} ---")
# Physical position: Coordinates from the center of the Sun (0, 0, 0)
print(f"Position (x, y, z) [km]: {position}")
# Distance from Sun in km and Astronomical Units (1 AU ≈ 149.6 million km)
print(f"Distance to Sun:        {dist_km:,.1f} km ({dist_km / 1.496e8:.3f} AU)")

# =============================================================================
# STEP 5: Query Orientation Kernel for Rotation Matrix & North Pole
# =============================================================================
# spice.pxform queries the planetary constants kernel (pck00011.tpc) to compute
# the 3x3 rotation matrix from inertial space (J2000) to the planet's rotating surface.
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

# =============================================================================
# STEP 6: Unload All Kernels from Memory (Cleanup)
# =============================================================================
spice.kclear()