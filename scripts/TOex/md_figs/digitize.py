# Digitise the carbonated-peridotite experiments of Dasgupta & Hirschmann (2010, EPSL, Fig. 5a):
# filled blue circles = subsolidus, open red circles = supersolidus. Output: P (GPa), T (degC), kind.
import numpy as np, matplotlib.pyplot as plt
from scipy import ndimage as nd
im = plt.imread('fig-000.png')[..., :3].astype(float)
H, W = im.shape[:2]
r, g, b = im[..., 0], im[..., 1], im[..., 2]
dark = (r + g + b) / 3 < 0.25
# panel (a) frame: the longest dark vertical and horizontal runs in the left half
colsum = dark[:, :W // 2].sum(axis=0); rowsum = dark[:, :W // 2].sum(axis=1)
xs = np.where(colsum > 0.6 * colsum.max())[0]; ys = np.where(rowsum > 0.6 * rowsum.max())[0]
x0, x1 = xs.min(), xs.max(); y0, y1 = ys.min(), ys.max()
print('frame px', x0, x1, y0, y1)
T = lambda x: 900 + (x - x0) / (x1 - x0) * 1000
P = lambda y: 0 + (y - y0) / (y1 - y0) * 35
inside = np.zeros_like(dark); inside[y0 + 3:y1 - 2, x0 + 3:x1 - 2] = True
mean = (r + g + b) / 3
blue = (b - r > 0.04) & (b - g > 0.04) & (mean < 0.6) & inside
red = (r - g > 0.12) & (r - b > 0.10) & inside
out = []
def blobs(mask, rmin):
    """centres of disk-like parts of a mask: points at least rmin px inside it (lines and text strokes are thinner)"""
    dist = nd.distance_transform_edt(mask)
    lab, n = nd.label(dist >= rmin)
    return nd.center_of_mass(dist >= rmin, lab, range(1, n + 1))
for cy, cx in blobs(blue, 3.0): out.append((P(cy), T(cx), 'sub'))
for cy, cx in blobs(nd.binary_fill_holes(nd.binary_closing(red, iterations=1)), 2.5): out.append((P(cy), T(cx), 'super'))
out.sort()
with open('dh2010_fig5a_digitised.csv', 'w') as f:
    f.write('# Dasgupta & Hirschmann (2010) EPSL 298, Fig. 5a, digitised from the published raster; P (GPa), T (degC)\n')
    f.write('# sub = subsolidus run (filled blue), super = supersolidus run (open red)\nP_GPa,T_C,kind\n')
    for p, t, k in out: f.write(f'{p:.2f},{t:.0f},{k}\n')
print(len(out), 'points:', sum(k == 'sub' for _, _, k in out), 'sub', sum(k == 'super' for _, _, k in out), 'super')
# check plot over the original
fig, ax = plt.subplots(figsize=(7, 6), dpi=110)
ax.imshow(im[:, :W // 2 + 20]);
for p, t, k in out:
    x = x0 + (t - 900) / 1000 * (x1 - x0); y = y0 + p / 35 * (y1 - y0)
    ax.plot(x, y, 'x', color='lime' if k == 'sub' else 'k', ms=6, mew=1.2)
plt.savefig('digitise_check.png')
# zoom on the crowded low-pressure corner
fig, ax = plt.subplots(figsize=(7, 5), dpi=110)
ax.imshow(im[77:330, 79:400])
for p, t, k in out:
    x = x0 + (t - 900) / 1000 * (x1 - x0) - 79; y = y0 + p / 35 * (y1 - y0) - 77
    if y < 253 and x < 321: ax.plot(x, y, '+', color='lime' if k == 'sub' else 'k', ms=7, mew=1)
plt.savefig('digitise_zoom.png')
