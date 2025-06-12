import matplotlib.pyplot as plt
import numpy as np

import postprocess as pp

steps = pp.available_steps("rho")
if not steps:
    raise SystemExit("no rho")
xs, ys = pp.read_grid("rho")
is1d = ys is None
nx = len(xs)
for s in steps:
    rho = pp.load_field("rho", s, xs, ys)
    plt.figure()
    if is1d:
        plt.plot(xs, rho)
        plt.xlabel('x'); plt.ylabel('rho')
    else:
        X, Y = np.meshgrid(xs, ys)
        plt.contourf(X, Y, rho, levels=40, cmap='viridis')
        plt.gca().set_aspect('equal')
    plt.title(f"rho step {s}")
    out = f"Result/rho_{s}.png"
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out}")
