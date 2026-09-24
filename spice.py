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

# Position & velocity relative to Sun (J2000 frame)
state, lt = spice.spkezr(target, et, "J2000", "NONE", "SUN")
position = state[:3]
dist_km = np.linalg.norm(position)

print(f"\n--- {planet.capitalize()} on {date_str} ---")
print(f"Position (x, y, z) [km]: {position}")
print(f"Distance to Sun:        {dist_km:,.1f} km ({dist_km / 1.496e8:.3f} AU)")

# Orientation (Transformation from J2000 to Body-Fixed frame)
rot_matrix = spice.pxform("J2000", frame, et)
pole = rot_matrix.T @ np.array([0.0, 0.0, 1.0])
_, ra, dec = spice.recrad(pole)

print(f"\nOrientation Frame:      {frame}")
print(f"North Pole in J2000:    RA = {np.degrees(ra):.2f}°, Dec = {np.degrees(dec):.2f}°")
print(f"Rotation Matrix (J2000 -> {frame}):")
print(f"  X = [{rot_matrix[0, 0]:.6f}, {rot_matrix[0, 1]:.6f}, {rot_matrix[0, 2]:.6f}]")
print(f"  Y = [{rot_matrix[1, 0]:.6f}, {rot_matrix[1, 1]:.6f}, {rot_matrix[1, 2]:.6f}]")
print(f"  Z = [{rot_matrix[2, 0]:.6f}, {rot_matrix[2, 1]:.6f}, {rot_matrix[2, 2]:.6f}]")

spice.kclear()