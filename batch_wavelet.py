import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from matplotlib.colors import LogNorm
from astropy.visualization import ImageNormalize, LogStretch

# cofigure paths
FITS_PATH     = '/Users/alanisalvarado/Documents/uprm/GRAD-MAP/galaxies/ngc2997/ngc2997_miri_lv3_f770w_i2d_anchor.fits'
#CUTOUTS_CSV   = '/Users/alanisalvarado/Documents/uprm/GRAD-MAP/galaxies/ngc4951/ngc4951_cutouts.csv'      
PHANGS_CSV    = '/Users/alanisalvarado/Documents/uprm/GRAD-MAP/phangs_doc.csv'           
GALAXY_NAME   = 'NGC2997'                       
OUTPUT_DIR    = 'fractal_results'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR,'plots'), exist_ok=True)

# creating individual folder
gal_dir   = os.path.join(OUTPUT_DIR, GALAXY_NAME)
gal_plots = os.path.join(gal_dir, 'plots')
os.makedirs(gal_plots, exist_ok=True)

#  load the FITS, find first 2D HDU
hdul = fits.open(FITS_PATH)
data = None
for h in hdul:
    if h.data is not None and getattr(h.data, 'ndim',0) >=2:
        data = np.squeeze(h.data)
        header = h.header
        break
if data is None:
    raise RuntimeError("No 2D image found in FITS.")
wcs = WCS(header)
pixscale = abs(header.get('CDELT1',1.0)) * 3600.0  

print(f"Loaded FITS: shape={data.shape}, pixscale={pixscale:.3f}\"/pix")

# read through spreadsheet
phangs = pd.read_csv(PHANGS_CSV)
row = phangs.loc[phangs['Name']==GALAXY_NAME]
if row.empty:
    raise RuntimeError(f"{GALAXY_NAME} not in PHANGS CSV.")
dist_mpc = float(row['Distance'].values[0])
dist_pc  = dist_mpc * 1e6
arcsec2pc = dist_pc / 206265.0
print(f"Distance = {dist_mpc:.2f} Mpc -> {arcsec2pc:.2f} pc/arcsec")

# read cutouts 
lines = []
# with open(CUTOUTS_CSV) as f:
#     for L in f:
#         L=L.strip()
#         if not L or L.startswith('#'): 
#                 continue
#    # pick up any line that contains centerbox or rotbox
#         if ('centerbox' in L.lower()) or ('rotbox' in L.lower()):
#                 lines.append(L)
# print(f"Found {len(lines)} region lines in {CUTOUTS_CSV}")

# helper to parse one line
def parse_line(line):
    # centerbox [[ra, dec], [sx, sy]]
    # rotbox    [[ra, dec], [sx, sy], angle]
    nums = re.findall(r'\[+\s*([0-9.+-]+)deg[, ]+([0-9.+-]+)deg\]\D*\[\s*([0-9.+-]+)arcsec\D*([0-9.+-]+)arcsec', line)
    if not nums:
        return None
    ra, dec, sx, sy = map(float, nums[0])
    # optional rotation
    rot = 0.0
    mrot = re.search(r',\s*([0-9.+-]+)deg\]$', line)
    if mrot: rot = float(mrot.group(1))
    return dict(ra=ra, dec=dec, sx=sx, sy=sy, rot=rot)

# wavelet power & lmax detection
def wavelet_power(img):
    h,w = img.shape
    max_scale = min(h,w)//4
    scales = 2**np.arange(0, np.log2(max_scale))
    scales = scales[scales>=1].astype(int)
    power = []
    for s in scales:
        x = np.linspace(-4*s,4*s,8*s)
        y = x[:,None]
        wave = np.exp(-(x**2+y**2)/(2*s*s))*(1 - (x**2+y**2)/(2*s*s))
        wave /= np.sqrt((wave**2).sum())
        conv = np.abs(np.fft.ifft2(np.fft.fft2(img)*np.fft.fft2(wave,img.shape)))
        power.append(np.mean(conv**2))
    return np.array(scales), np.array(power)

def detect_lmax(scales, power, drop=0.5):
    r = power[1:]/power[:-1]
    idx = np.where(r<drop)[0]
    return scales[idx[0]] if idx.size else scales.max()

# loop regions
# results = []
# for i,line in enumerate(lines, start=1):
#     reg = parse_line(line)
#     if reg is None:
#         print(f"  skip unparsable: {line}")
#         continue
#     print(f"\n-- region {i}: RA={reg['ra']:.4f}, Dec={reg['dec']:.4f}, size={reg['sx']}x{reg['sy']}\" --")
#     # world->pix
#     xpix, ypix = wcs.world_to_pixel_values(reg['ra'], reg['dec'])
#     dx = int(reg['sx']/pixscale/2)
#     dy = int(reg['sy']/pixscale/2)
#     x0,x1 = max(0,int(xpix-dx)), min(data.shape[1], int(xpix+dx))
#     y0,y1 = max(0,int(ypix-dy)), min(data.shape[0], int(ypix+dy))
#     cut = data[y0:y1, x0:x1].copy()
#     cut[np.isnan(cut)] = 0
#     cut[cut<0] = 0
#     print(f"  cutout shape={cut.shape}")

#     if cut.size < 16:
#         print("  skip: too small")
#         continue

#     sc, pw = wavelet_power(cut)
#     lmax_pix = detect_lmax(sc,pw)
#     mask = (sc>=1) & (sc<=lmax_pix)
#     sc_fit = sc[mask]*pixscale*arcsec2pc
#     pw_fit = pw[mask]
#     if len(sc_fit)<3 or np.any(pw_fit<=0):
#         print("  skip fit: insufficient")
#         slope, intercept, r2 = np.nan, np.nan, np.nan
#     else:
#         ls = np.log10(sc_fit)
#         lp = np.log10(pw_fit)
#         slope, intercept = np.polyfit(ls,lp,1)
#         # R^2
#         pred = slope*ls+intercept
#         ssr = np.sum((lp-pred)**2)
#         sst = np.sum((lp-lp.mean())**2)
#         r2 = 1-ssr/sst

#     # plot
#     # plot: show cutout + power spectrum side by side
#     fig, (ax_img, ax_pow) = plt.subplots(1, 2, figsize=(12, 5))

#     # left panel: the cutout image
#     norm = ImageNormalize(vmin=cut.min(), vmax=cut.max(),
#                         stretch=LogStretch()) 
#     im = ax_img.imshow(cut, origin='lower', cmap='viridis', norm=norm)
#     ax_img.set_title(f"Region {i:02d}")
#     ax_img.axis('off')
#     fig.colorbar(im, ax=ax_img, fraction=0.046, pad=0.04, label='Intensity')

#     # right panel: wavelet power spectrum
#     sc_pc = sc * pixscale * arcsec2pc
#     ax_pow.loglog(sc_pc, pw, 'o-', label='Wavelet Power')
#     if not np.isnan(slope):
#         fit_line = 10**(slope * np.log10(sc_fit) + intercept)
#         ax_pow.loglog(sc_fit, fit_line, 'r--', label=f'slope={slope:.2f}')
#     ax_pow.axvline(lmax_pix * pixscale * arcsec2pc, ls='--', color='gray',
#                    label=f'lmax={lmax_pix * pixscale * arcsec2pc:.1f} pc')
#     ax_pow.set_xlabel('Scale (pc)')
#     ax_pow.set_ylabel('Wavelet Power')
#     ax_pow.set_title('Power Spectrum')
#     ax_pow.grid(which='both', linestyle='--', alpha=0.7)
#     ax_pow.legend()

#     # save figure
#     out = os.path.join(gal_plots, f'region_{i:02d}.png')
#     fig.savefig(out, dpi=200, bbox_inches='tight')
#     plt.close(fig)
#     print(f"  saved {out}")


#     results.append({
#         'region':i,
#         'slope':slope,
#         'r2':r2,
#         'lmax_pc':lmax_pix*pixscale*arcsec2pc
#     })

    # full image: show whole galaxy + its power spectrum
print("\n-- full image --")
cut_all = np.clip(np.nan_to_num(data), 0, None)
sc_all, pw_all = wavelet_power(cut_all)
lmax_all = detect_lmax(sc_all, pw_all)
mask_all = (sc_all >= 1) & (sc_all <= lmax_all)
sc_fit_all = sc_all[mask_all] * pixscale * arcsec2pc
pw_fit_all = pw_all[mask_all]
slope_all, intercept_all = np.nan, np.nan
if len(sc_fit_all) >= 3 and not np.any(pw_fit_all <= 0):
    ls_all = np.log10(sc_fit_all)
    lp_all = np.log10(pw_fit_all)
    slope_all, intercept_all = np.polyfit(ls_all, lp_all, 1)

# two‐panel plot for full image
fig, (ax_img, ax_pow) = plt.subplots(1, 2, figsize=(12, 5))
# left: full galaxy
norm = ImageNormalize(vmin=0.1, vmax=50,
                    stretch=LogStretch()) 
im = ax_img.imshow(cut_all, origin='lower', cmap='viridis', norm=norm)
ax_img.set_title('NGC2997')
ax_img.axis('off')
fig.colorbar(im, ax=ax_img, fraction=0.046, pad=0.04, label='Intensity')

# right: its power spectrum
sc_pc_all = sc_all * pixscale * arcsec2pc
ax_pow.loglog(sc_pc_all, pw_all, 'o-', label='Wavelet Power')
if not np.isnan(slope_all):
    fit_line_all = 10**(slope_all * np.log10(sc_fit_all) + intercept_all)
    ax_pow.loglog(sc_fit_all, fit_line_all, 'r--', label=f'slope={slope_all:.2f}')
ax_pow.axvline(lmax_all * pixscale * arcsec2pc, ls='--', color='gray',
                label=f'lmax={lmax_all * pixscale * arcsec2pc:.1f} pc')
ax_pow.set_xlabel('Scale (pc)')
ax_pow.set_ylabel('Wavelet Power')
ax_pow.set_title('Full Image Spectrum')
ax_pow.grid(which='both', linestyle='--', alpha=0.7)
ax_pow.legend()

# save full-image figure
out_all = os.path.join(gal_plots, 'full_image.png')
fig.savefig(out_all, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f"  saved {out_all}")

# write summary CSV
# pd.DataFrame(results).to_csv(os.path.join(OUTPUT_DIR,'results.csv'),index=False)
# print("\nDone.  Results in", OUTPUT_DIR)
    # write full-galaxy spectrum CSV with exact (scale, power) points
sc_pc_all = sc_all * pixscale * arcsec2pc
df_points = pd.DataFrame({
    'scale_pc': sc_pc_all,
    'power': pw_all
})
out_csv = os.path.join(gal_dir, f'{GALAXY_NAME}_spectrum.csv')
df_points.to_csv(out_csv, index=False)
print(f"Saved full-galaxy spectrum points CSV: {out_csv}")



