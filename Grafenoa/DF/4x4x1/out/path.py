import numpy as np

points = [
    np.array([0.000000, 0.000000, 0.0]),  # Gamma
    np.array([0.500000, 0.000000, 0.0]),  # M
    np.array([0.333333, 0.333333, 0.0]),  # K
    np.array([0.000000, 0.000000, 0.0]),  # Gamma
]
npts = 50  # igual que en matdyn.in

path = []
for i in range(len(points) - 1):
    segment = np.linspace(points[i], points[i+1], npts, endpoint=False)
    path.extend(segment)
path.append(points[-1])

np.savetxt('path.dat', np.array(path))
print(f"path.dat generado: {len(path)} puntos")  # debe dar 151
