#SLM_mode_decomp
# Imports
from HEDS import *
from HEDS import holoeye_slmdisplaysdk_slm as heds
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hermite, factorial
from scipy.io import loadmat
import math
import time
from pypylon import pylon
import cv2
import csv
from scipy.optimize import curve_fit

def gaussian(x, A, mu, sigma):
    return A * np.exp(- (x - mu)**2 / (2 * sigma**2))


# Connecting to SLM and Camera

def connect_slm():
    heds.SDK.Init(4,1)
    slm = heds.SLM.Init(openPreview=True, previewScale=1)

    # Check for errors
    if slm.errorCode() != heds.HEDSERR_NoError:
        print(f"Error connecting to SLM: {heds.SDK.errorString(slm.errorCode())}")
    else:
        # 3. If the connection is successful, you can now display data on the SLM.
        #    Here, we'll just show a blank screen.
        slm.showBlankScreen(000) # Show a mid-gray screen
        print("connected")
    return slm

def connect_camera():
    tl_factory = pylon.TlFactory.GetInstance()
    devices = tl_factory.EnumerateDevices()
    if not devices:
        raise RuntimeError("No Basler camera detected. Check connection and drivers.")
    
    camera = pylon.InstantCamera(tl_factory.CreateFirstDevice())
    camera.Open()  # Open the camera immediately
    print("Connected to:", camera.GetDeviceInfo().GetModelName())
    return camera

# functions to capture and process camera images

def capture_single_image(camera, exposure_time):
    if not camera.IsOpen():
        camera.Open()
    camera.PixelFormat.SetValue("Mono8")
    camera.ExposureTime.SetValue(exposure_time)

    camera.StartGrabbingMax(1)
    grabResult = camera.RetrieveResult(5000)

    if grabResult.GrabSucceeded():
        img = grabResult.Array  # ✅ This is a NumPy array
    else:
        grabResult.Release()
        raise RuntimeError("Image capture failed.")

    grabResult.Release()
    return img


def subtract_images(img1, img2):
    # Ensure same shape
    if img1.shape != img2.shape:
        raise ValueError("Images must have the same shape to subtract.")
    diff = np.abs(img2.astype(np.int16) - img1.astype(np.int16))
    # Clip to avoid negative values and convert back to uint8
    diff = np.clip(diff, 0, 255).astype(np.uint8)
    # Alternatively, if you want to keep only positive differences:
    # Uncomment the next line to use it instead of the above line
    #diff = np.clip(img2.astype(np.int16) - img1.astype(np.int16), 0, None).astype(np.uint8)
    return diff

def intensity_calculation(image, type = 'mean'):
    # Calculate the intensity of the image
    if type == 'mean':
        intensity = np.mean(image)

    if type == 'sum':
        intensity = np.sum(image)
    return intensity


# knife_edge and mode decomp grating patterns:

import numpy as np

import numpy as np

def partial_grating_phase(
    resolution=(1920, 1080),
    grating_period=100,
    coverage_fraction=0.5,
    direction='horizontal',
    constant_phase=0.0,
):
    """
    Create a phase-only diffraction grating that covers only part of the SLM.

    Parameters:
    - resolution: (width, height) of SLM in pixels
    - grating_period: pixels per 2π phase cycle
    - coverage_fraction: portion of SLM to cover with grating (0 to 1)
    - direction: 'horizontal' (vertical grating lines, vary along y, coverage top→bottom)
                 'vertical'   (horizontal grating lines, vary along x, coverage left→right)
    - constant_phase: phase value outside the grating region
    """

    width, height = resolution
    phase_pattern = np.ones((height, width)) * constant_phase

    if direction == 'horizontal':
        # Vertical grating lines → vary along y
        y = np.arange(height)
        grating = 2 * np.pi * (y % grating_period) / grating_period
        grating = np.tile(grating[:, np.newaxis], (1, width))

        # Coverage grows top → bottom
        cover_height = int(height * coverage_fraction)
        phase_pattern[:cover_height, :] = grating[:cover_height, :]

    elif direction == 'vertical':
        # Horizontal grating lines → vary along x
        x = np.arange(width)
        grating = 2 * np.pi * (x % grating_period) / grating_period
        grating = np.tile(grating, (height, 1))

        # Coverage grows left → right
        cover_width = int(width * coverage_fraction)
        phase_pattern[:, :cover_width] = grating[:, :cover_width]

    return phase_pattern





import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import hermite, j1  # <-- use J1 (first kind, order 1)

def generate_hg_hologram(m, n, w0=1.0, center=(0.0, 0.0), theta_deg=0,
                         conjugate=True, plot_output=True,
                         Grating=False, Grating_period=(100, 100)):
    """
    Phase-only hologram for HG(m,n) using the *first-order* Bessel mapping:
        A_target in [0,1]  ->  F in [0, x0]  s.t.  J1(F) = A_target * J1_max
    where x0 ≈ 1.84118, J1_max ≈ 0.581865.

    Returns:
        slm_mask : uint8 2D array in [0,255]
    """
    # --- SLM grid ---
    H, V = 1920, 1080
    pixel_size = 8e-3  # mm
    x_1d = np.arange(-H/2, H/2) * pixel_size
    y_1d = np.arange(-V/2, V/2) * pixel_size
    X, Y = np.meshgrid(x_1d, y_1d)

    # Center + rotate
    x_center, y_center = center
    Xu, Yu = X - x_center, Y - y_center
    th = np.deg2rad(theta_deg)
    X =  np.cos(th)*Xu - np.sin(th)*Yu
    Y =  np.sin(th)*Xu + np.cos(th)*Yu

    # --- Beam parameters ---
    lmbda = 633e-6  # mm
    z = 1e-5        # mm (near SLM plane)
    k = 2*np.pi/lmbda
    zr = np.pi*w0**2/lmbda
    w  = w0*np.sqrt(1 + (z/zr)**2)
    R  = z*(1 + (zr/z)**2) if z != 0 else np.inf

    # Hermite polynomials (physicists')
    Hm = hermite(m)
    Hn = hermite(n)
    Hx = Hm(np.sqrt(2)*X/w)
    Hy = Hn(np.sqrt(2)*Y/w)

    # Normalization, Gouy, curvature
    rc   = np.sqrt(2**(1-n-m)/(np.pi*math.factorial(n)*math.factorial(m))) / w
    gouy = (n + m + 1)*np.arctan(z/zr)
    rho2 = X**2 + Y**2

    HG = (rc*Hx*Hy *
          np.exp(1j*gouy) *
          np.exp(-rho2/w**2) *
          np.exp(-1j*k*rho2/(2*R)) *
          np.exp(1j*k*z))

    if conjugate:
        HG = np.conj(HG)

    # Normalize power
    HG = HG / np.sqrt(np.sum(np.abs(HG)**2))

    # Amplitude A ∈ [0,1], phase φ
    A  = np.abs(HG)
    A /= (A.max() if A.max() != 0 else 1.0)
    A  = np.clip(A, 0.0, 1.0)
    ph = np.angle(HG)

    # -------- Bessel J1 inversion (robust + fast, no root-finding) ----------
    X0      = 1.8411837813406593   # arg at max of J1
    J1_MAX  = 0.5818652247353657   # max J1 value
    # Build a high-res monotone table for J1 on [0, X0]
    _xs = np.linspace(0.0, X0, 4096)
    _ys = j1(_xs)
    # Numerically-enforce monotonic increase (safety against tiny ripples)
    _ys_mono = np.maximum.accumulate(_ys)

    # Map target amplitude to J1 range and invert by interpolation
    target = A * J1_MAX                         # now in [0, J1_MAX]
    F = np.interp(target, _ys_mono, _xs)        # F in [0, X0]
    # ------------------------------------------------------------------------

    # Carrier grating
    if Grating:
        nx, ny = Grating_period
    else:
        nx, ny = 0, 0
    H_mm, V_mm = H*pixel_size, V*pixel_size
    gx, gy = nx/H_mm, ny/V_mm
    Hol = F * np.sin(ph + 2*np.pi*(gx*X + gy*Y))

    # Normalize to 8-bit
    Hol -= Hol.min()
    slm_mask = (Hol / (Hol.max() if Hol.max() != 0 else 1.0) * 255).astype(np.uint8)

    if plot_output:
        plt.figure(figsize=(12, 6.75))
        plt.imshow(slm_mask, cmap='gray', vmin=0, vmax=255)
        plt.title(f'HG({m},{n}), w0={w0} mm, Conjugate={conjugate}', fontsize=18)
        plt.axis('off')
        plt.show()

    return slm_mask

# Performing the knife edge to find beam width and position
def digital_knife_edge(SLM, camera, orientation, exposure_time, steps = 500, i_j = 1, grating_period = 3):
    """plots a knife edge intensity pattern

    process:
    1. use partial grating to create a moving phase pattern on SLM
    2. take intensity measurements by image diff at each step
    3. plot the intensity profile
    """
    # Initialize SLM and camera
    #slm = heds.SLM.Init(openPreview=True, previewScale=0.5)
    # camera = pylon.InstantCamera(pylon.TlFactory.GetInstance().CreateFirstDevice())
    # camera.Open()
    # camera.PixelFormat.SetValue("Mono8")
    # camera.ExposureTime.SetValue(155000)
    # camera.close()

    x = np.linspace(0, 1, steps)
    intensity_profile = []
    img1 = capture_single_image(camera, exposure_time=exposure_time)  # Capture the initial image (flush)
    for dist in x:
        #print(dist)
        # Update SLM with partial grating phase
        SLM.showPhaseData(partial_grating_phase(
            resolution=(1920, 1080),
            grating_period= grating_period,
            coverage_fraction=dist,  # Varying coverage fraction
            direction= 'vertical' if orientation == 'vertical' else 'horizontal',
            #save_image=False
        ))
        #time.sleep(0.1)  # Allow SLM to update
        # Capture image
        img2 = capture_single_image(camera, exposure_time=exposure_time)
        # Subtract initial image to get intensity change
        diff = subtract_images(img1, img2)
        if i_j:
            img1 = img2  # Update img1 for the next iteration
        # Calculate intensity and store it
        intensity = intensity_calculation(diff)
        intensity_profile.append(intensity)
    #SLM.showBlankScreen(0)
    plt.plot(x, intensity_profile)
    plt.show()

    return intensity_profile

def find_beam_properties(slm, camera, exposure_time):
    """
    Find beam location and sigma_x, sigma_y using digital knife edge
    """
    # Perform knife edge scans in both directions

    timestamp = int(time.time())

    intensity_x = digital_knife_edge(slm, camera, orientation='horizontal', exposure_time=exposure_time)
    intensity_y = digital_knife_edge(slm, camera, orientation='vertical', exposure_time=exposure_time)

    # Save intensity profiles as CSV with timestamp in filename
    np.savetxt(f'intensity_x_{timestamp}.csv', np.array(intensity_x), delimiter=',')
    np.savetxt(f'intensity_y_{timestamp}.csv', np.array(intensity_y), delimiter=',')

    # processing the measured intensity profiles
    
    #Fitting gaussian in horizontal orientation
    x = np.linspace(1,1920,500)
    popt_x, _ = curve_fit(gaussian, x, intensity_x, p0=[np.max(intensity_x), 960, 100])
    sigma_x = popt_x[2]
    mu_x = popt_x[1]
    #Fitting gaussian in vertical orientation
    y = np.linspace(1,1080,500)
    popt_y, _ = curve_fit(gaussian, y, intensity_y, p0=[np.max(intensity_y), 960, 100])
    sigma_y = popt_y[2]
    mu_y = popt_y[1]

    # Optionally, you can visualize the fits
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.plot(x, intensity_x, label='Measured Intensity (X)')
    plt.plot(x, gaussian(x, *popt_x), label='Gaussian Fit (X)', linestyle='--')
    plt.title('Horizontal Intensity Profile')
    plt.xlabel('Pixel')
    plt.ylabel('Intensity')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(y, intensity_y, label='Measured Intensity (Y)')
    plt.plot(y, gaussian(y, *popt_y), label='Gaussian Fit (Y)', linestyle='--')
    plt.title('Vertical Intensity Profile')
    plt.xlabel('Pixel')
    plt.ylabel('Intensity')
    plt.legend()

    plt.tight_layout()
    plt.show()

    plt.savefig(f'beam_profiles_{timestamp}.png')

    return [mu_x, sigma_x, mu_y, sigma_y]
