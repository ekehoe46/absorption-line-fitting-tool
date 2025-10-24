import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from glob import glob
import argparse
import os
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

################# CONSTANTS (CHANGE AS NEEDED #################
SIM_NUM = 10 #number of simulated spectra you are creating
EDGE = 35 #interval width to look for abs line (in angstroms)
GUESS = 10 #the allowed range around the expected line center where the fit searches for the Gaussian peak
FLUX_MEAN_SIG = 4 #minimum absorption line significance level to be detected
FILE_NAME = 'combined_absorption_results.txt'
###############################################################

def extract_spectrum(filepath: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reads a text spectrum with columns: wavelength, flux, error, and normalizes the flux.

    Parameters:
        filepath (str): Path to the spectrum file.

    Returns:
        wave (np.ndarray): Wavelength array.
        flux_norm (np.ndarray): Normalized flux array.
        err_norm (np.ndarray): Normalized error array.
    """
    data = pd.read_csv(filepath, sep=r'\s+')
    wave = data['lam'].values
    flux = data['flux'].values
    err = data['sigma'].values
    norm = np.nanmean(flux)

    return wave, flux / norm, err / norm


def get_spectrum_files(directory: str) -> list[str]:
    """
    Returns a sorted list of all .txt spectrum files in the specified directory.

    Parameters:
        directory (str): Path to the folder containing spectrum files.

    Returns:
        List[str]: Sorted list of full paths to .txt files.
    """
    return sorted(glob(os.path.join(directory, "*.txt")))


def generate_gaussian_realizations(flux: np.ndarray, flux_err: np.ndarray, 
                                   sim_num: int = SIM_NUM) -> np.ndarray:
    """
    Generates Gaussian realizations of a flux array based on its errors
    for Monte Carlo error propagation.

    Parameters:
        flux (np.ndarray): Original flux array.
        flux_err (np.ndarray): 1D array of flux uncertainties.
        sim_num (int): Number of simulated spectra to generate.

    Returns:
        np.ndarray: Array of shape (sim_num, len(flux)) with simulated fluxes.
    """
    return np.array([
        [flux[i] + np.random.normal(0, flux_err[i]) for i in range(len(flux))]
        for _ in range(sim_num)
    ])

def fit_absorption_continuum(blue1: float, blue2: float, red1: float, red2: float, 
                             wave: np.ndarray, flux: np.ndarray) -> tuple[float, float]:
    """
    Estimates the local continuum level around an absorption line
    by averaging flux in blue and red wavelength windows.

    Parameters:
        blue1/red1, blue2/red2 (float): Lower and upper limits of the blue/red continuum window (in angstroms).
        wave (np.ndarray): Wavelength array
        flux (np.ndarray): Flux array

    Returns:
        tuple[float, float]:
            slope (float): Continuum slope (set to 0 for a flat continuum).
            intercept (float): Mean continuum level from blue and red windows.
    """
    #Calculate average continuum from blue and red line sections
    c_blue = np.mean(flux[(wave >= blue1) & (wave <= blue2)])
    c_red = np.mean(flux[(wave >= red1) & (wave <= red2)])
    #assume slope is 0 to have a perfectly flat line
    return 0, 0.5 * (c_red + c_blue) 

def integrate_flux_line(wave: np.ndarray, flux: np.ndarray, 
                        ion_wave: float, 
                        blue1: float, blue2: float, red1: float, red2: float, 
                        z: float) -> float:
    """
    Integrates the continuum subtracted flux around a target absorption line.

    Parameters:
        wave (np.ndarray): Wavelength array
        flux (np.ndarray): Flux array
        ion_wave (float): Rest wavelength of the target absorption line (in angstroms).
        blue1/red1, blue2/red2 (float): Blue/red continuum window limits (in angstroms).
        z (float): Redshift of the object.

    Returns:
        float: Integrated flux of the target line (negative for absorption). Returns NaN if invalid.
    """

    delta_lambda = np.mean(np.diff(wave))
    ion_wave_z = ion_wave * (1 + z)
    lower_bound = ion_wave_z - EDGE
    upper_bound = ion_wave_z + EDGE

    #skips lines not in spectrum
    if lower_bound < wave.min() or upper_bound > wave.max():
        return np.nan 

    mask = (wave >= lower_bound) & (wave <= upper_bound)
    if np.sum(mask) == 0:
        return np.nan

    wave_cut, flux_cut = wave[mask], flux[mask]

    #finds continuum if region exists in spectrum
    try:
        m, b = fit_absorption_continuum(
            blue1 * (1 + z), blue2 * (1 + z),
            red1 * (1 + z), red2 * (1 + z),
            wave, flux
        )
    except:
        return np.nan

    #calculating integrated flux
    continuum = m * wave_cut + b
    residual = flux_cut - continuum
    integrated_flux = np.sum(residual) * delta_lambda

    return integrated_flux

#finds area under the curve for where the line should be to determine if the line is significantlly in the spectrum
def non_parametric_fit(wave: np.ndarray, spec_realizations: np.ndarray, 
                       line_file: pd.DataFrame, z: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculates the non-parametric fluxes for absorption lines using area under the curve for each realization to determine
    if the line is signficantly detected.

    Parameters:
        wave (np.ndarray): Wavelength array
        spec_realizations (np.ndarray): 2D array of shape (num_realizations, len(wave))
                                       with simulated flux spectra
        line_file (pd.DataFrame): DataFrame containing absorption line info with columns:
                                  'lam0', 'blue-window1', 'blue-window2',
                                  'red-window1', 'red-window2'.
        z (float): Redshift of the object

    Returns:
        tuple[np.ndarray, np.ndarray]:
            flux_mean (np.ndarray): Mean integrated flux for each line across all realizations
            flux_std (np.ndarray): Standard deviation of integrated flux for each line
    """ 

    #list to hold fluxes for each realization
    flux_arr = []
    
    #iterates over each simulated spectrum
    for sim_flux in spec_realizations:
        sim_flux_list= []

        #iterates over each line
        for idx, ion in enumerate(line_file['lam0']):
            blue1, blue2 = line_file['blue-window1'][idx], line_file['blue-window2'][idx]
            red1, red2 = line_file['red-window1'][idx], line_file['red-window2'][idx]

            flux_val = integrate_flux_line(wave, sim_flux, ion, blue1, blue2, red1, red2, z)
            sim_flux_list.append(flux_val) 

        flux_arr.append(sim_flux_list) 

    flux_arr = np.array(flux_arr)
    #mean fluxes found for each line by averaging fluxes for each line with associated standard deviations 
    flux_mean = np.nanmean(flux_arr, axis=0) 
    flux_std = np.nanstd(flux_arr, axis=0) 

    return flux_mean, flux_std

def absorption_model(wave: np.ndarray,
                      intercept: float, slope: float, 
                      amplitude: float, width: float, centroid: float) -> np.ndarray:
    """
    Gaussian absorption line model with a linear continuum. 
    The functional form used to model the absorption lines.

    Parameters:
        wave(np.ndarray): Wavelength array where the model is evalulated
        intercept (float): Linear continuum offset (baseline flux)
        slope (float): Linear continuum slope (usually 0)
        amplitude (float): Amplitude of the Gaussian absorption (negative for absorption lines)
        width (float): Standard deviation of the Gaussian profile (in angstroms)
        centroid (float): Center of the Gaussian line (observed-frame wavelength in angstroms)

    Returns:
        np.ndarray: Model flux values at each wavelength
    """
    gaussian = amplitude * np.exp(-0.5 * ((wave - centroid) / width) ** 2)
    
    return intercept + slope * wave + gaussian

def fit_line_curve(wave: np.ndarray, flux: np.ndarray, 
                   ion_wave_z: float, 
                   intercept: float, guess: float=GUESS) -> tuple[float, list[float]]:
    """
    Fits a Gaussian absorption line model to a section of the spectrum centered around a target absorption line.

    Parameters:
        wave(np.ndarray): Wavelength array in angstroms
        flux (np.ndarray): Flux array
        ion_wave_z (float): Wavelength of the target absorption line in the observed frame of the galaxy
        intercept (float): Inital guess for continuum level based on the average continuum flux on the blue and red sides of the target absorption line
        guess (float): Allowed wavelength range for the fitted centroid of the absorption line

    Returns:
         tuple:
            fitted_flux (float): Integrated Gaussian absorption line flux 
                                 Returns NaN if fit fails
            params (list[float]): Best fit parameters [intercept, slope, amplitude, width, centroid]
                                 Returns list of NaN values if fit fails

    """
    lower_bound = ion_wave_z - EDGE
    upper_bound = ion_wave_z + EDGE

    if lower_bound < wave.min() or upper_bound > wave.max():
        return np.nan, [np.nan]*5

    mask = (wave >= lower_bound) & (wave <= upper_bound)
    if np.sum(mask) == 0:
        return np.nan, [np.nan]*5
    
    idx_centroid = np.abs(wave - ion_wave_z).argmin()
    flux_at_centroid = flux[idx_centroid]
    
    # If flux is 0, skip fitting
    if np.abs(flux_at_centroid) == 0.0:
        return np.nan, [np.nan]*5

    wave_cut, flux_cut = wave[mask], flux[mask]

    try:
        #edit the params as necessary
        popt, _ = curve_fit(
            absorption_model, wave_cut, flux_cut,
            p0=(intercept, 0., -0.1, 1.5, ion_wave_z),
            bounds=([intercept, 0, -20., 0.25, ion_wave_z - guess],
                    [15., 1e-11, -0.01, 20., ion_wave_z + guess])
        )
        flux_fit = popt[2] * popt[3] * np.sqrt(2 * np.pi)
        return flux_fit, popt
    except:
        return np.nan, [np.nan]*5


def curvefit_measurements(wave: np.ndarray, flux_realizations: np.ndarray, 
                          line_file: pd.DataFrame, 
                          z: float, 
                          flux_mean: np.ndarray, flux_std: np.ndarray) -> tuple[list, ...]:
    """
    Runs Gaussian line fitting on each simulated spectrum to estimate mean and uncertainty for each fit parameter for detected absorption lines.

    A line is considered detected when:
        abs(flux_mean) > significance_level * flux_std
        and flux_mean < 0  (indiciative of absorption)

    Parameters:
        wave (np.ndarray): Wavelength array
        flux_realizations (np.ndarray): Array of simulated spectra
        line_file (pd.DataFrame): Table with columns defining lines and continuum windows.
        z (float): Redshift of the object
        flux_mean (np.ndarray): Mean non-parametric detection flux
        flux_std (np.ndarray): Standard deviation for detection flux

    Returns:
        tuple of lists. Each list has one entry per absorption line:
            flux_fit_mean, flux_fit_std,
            a_mean, a_std,
            b_mean, b_std,
            d_mean, d_std,
            s_mean, s_std,
            mu_mean, mu_std
    """   
    flux_fit_mean, flux_fit_std = [], []
    a_mean, a_std, b_mean, b_std, d_mean, d_std, s_mean, s_std, mu_mean, mu_std = [], [], [], [], [], [], [], [], [], []

    for idx, ion in enumerate(line_file['lam0']):
        if abs(flux_mean[idx]) > FLUX_MEAN_SIG * flux_std[idx] and flux_mean[idx] < 0:
            flux_vals = []
            a_list, b_list, d_list, s_list, mu_list = [], [], [], [], []

            for sim_flux in flux_realizations:
                m, b_c = fit_absorption_continuum(
                    line_file['blue-window1'][idx]*(1+z), line_file['blue-window2'][idx]*(1+z),
                    line_file['red-window1'][idx]*(1+z), line_file['red-window2'][idx]*(1+z),
                    wave, sim_flux
                )
                f, popt = fit_line_curve(wave, sim_flux, ion*(1+z), b_c)
                flux_vals.append(f)
                a_list.append(popt[0])
                b_list.append(popt[1])
                d_list.append(popt[2])
                s_list.append(popt[3])
                mu_list.append(popt[4])

            flux_fit_mean.append(np.nanmean(flux_vals))
            flux_fit_std.append(np.nanstd(flux_vals))

            a_mean.append(np.nanmean(a_list)); a_std.append(np.nanstd(a_list))
            b_mean.append(np.nanmean(b_list)); b_std.append(np.nanstd(b_list))
            d_mean.append(np.nanmean(d_list)); d_std.append(np.nanstd(d_list))
            s_mean.append(np.nanmean(s_list)); s_std.append(np.nanstd(s_list))
            mu_mean.append(np.nanmean(mu_list)); mu_std.append(np.nanstd(mu_list))
        else:
            flux_fit_mean.append(np.NaN); flux_fit_std.append(np.NaN)
            a_mean.append(np.NaN); a_std.append(np.NaN)
            b_mean.append(np.NaN); b_std.append(np.NaN)
            d_mean.append(np.NaN); d_std.append(np.NaN)
            s_mean.append(np.NaN); s_std.append(np.NaN)
            mu_mean.append(np.NaN); mu_std.append(np.NaN)

    return (flux_fit_mean, flux_fit_std,
            a_mean, a_std, b_mean, b_std, d_mean, d_std, s_mean, s_std, mu_mean, mu_std)


def process_spectrum(filepath: str, 
                     line_file: pd.DataFrame, redshift_df: pd.DataFrame) -> dict[str, any]:
    """
    Process a single spectrum by generating simulated realizations, 
    then performing a non parametric integration to find signficantly detected lines,
      and lastly running Gaussian line fits on the significantly detected lines.

    Parameters:
        filepath (str): Path to spectrum text file
        line_file (pd.DataFrame): Table defining absorption features and continuum estimation windows
        redshift_df (pd.DataFrame): Table with object names and redshift

    Returns:
        dict: Dictionary of lists with measurement results for each line
    """
    results = {
        'flux_mean': [], 'flux_std': [],
        'flux_fit_mean': [], 'flux_fit_std': [],
        'a_mean': [], 'a_std': [], 'b_mean': [], 'b_std': [],
        'd_mean': [], 'd_std': [], 's_mean': [], 's_std': [],
        'mu_mean': [], 'mu_std': []
    }
    obj = os.path.basename(filepath).replace('.txt', '')
    z = float(redshift_df.loc[redshift_df['Object'] == obj, 'z'].values[0])

    try:
        wave, base_flux, base_err = extract_spectrum(filepath)
        flux_realizations = generate_gaussian_realizations(base_flux, base_err)
        flux_mean, flux_std = non_parametric_fit(wave, flux_realizations, line_file, z)
        results['flux_mean'] = flux_mean
        results['flux_std'] = flux_std

        (f_fit_mean, f_fit_std, 
         a_mean, a_std, b_mean, b_std, d_mean, d_std, s_mean, s_std, mu_mean, mu_std) = curvefit_measurements(
             wave, flux_realizations, line_file, z, flux_mean, flux_std
         )
        results['flux_fit_mean'] = f_fit_mean; results['flux_fit_std'] = f_fit_std
        results['a_mean'] = a_mean; results['a_std'] = a_std
        results['b_mean'] = b_mean; results['b_std'] = b_std
        results['d_mean'] = d_mean; results['d_std'] = d_std
        results['s_mean'] = s_mean; results['s_std'] = s_std
        results['mu_mean'] = mu_mean; results['mu_std'] = mu_std

    except Exception as e:
        print(f"Error processing ID {obj}: {e}")
        n_lines = len(line_file)
        for key in results.keys():
            results[key] = [np.nan] * n_lines

    return results

def get_filename(output_dir: str, base_name: str = FILE_NAME) -> str:
    """
    Generates a versioned output filename that avoids overwriting existing files.

    Parameters:
        output_dir (str): Folder where the output should be saved
        base_name (str): Desired base filename

    Returns:
        str: A filepath inside output_dir with a version suffix if needed.
    """
    if not os.path.exists(base_name):
        return base_name

    base, ext = os.path.splitext(base_name)
    version = 1
    while True:
        new_name = f"./{output_dir}/{base}_v{version}{ext}"
        if not os.path.exists(new_name):
            return new_name
        version += 1

def main(spectra_dir: str, 
         line_file: str, 
         redshift_file: str, 
         output_dir: str) -> None:
    """
    Process all spectra in a directory and measures desired absorption lines

    Parameters:
        spectra_dir (str): Directory containing spectrum files (.txt)
        line_file (str): Path to absorption line list file
        redshift_file (str): Path to redshift table
        output_dir (str): Output folder for results

    Returns:
        No outputs, but a results file is written to the output folder.
    """

    os.makedirs(output_dir, exist_ok=True)
    redshift_df = pd.read_csv(redshift_file, sep=r'\s+', dtype={'Object': str})
    line_df = pd.read_csv(line_file, sep=r'\s+')
    spectra_files = get_spectrum_files(spectra_dir)
    if not spectra_files:
        print("No .txt spectra found in the provided directory.")
        return

    all_data = []
    for spec_file in spectra_files:
        print(f"Fitting lines to: {os.path.basename(spec_file)}")
        df = process_spectrum(spec_file, line_df, redshift_df)
        if df is not None:
            all_data.append(df)

    data = {
        'Object': [], 'z': [], 'Line': [],
        'NonParametric_Flux': [], 'NonParametric_Flux_std': [],
        'Gaussian_Flux': [], 'Gaussian_Flux_std': [],
        'a': [], 'a_std': [], 'b': [], 'b_std': [],
        'd': [], 'd_std': [], 's': [], 's_std': [], 'mu': [], 'mu_std': []
    }
    
    for i, obj in enumerate(redshift_df['Object']):
        for j, line in enumerate(line_df['linename']):
            res = all_data[i]
            data['Object'].append(obj)
            data['z'].append(redshift_df['z'][i])
            data['Line'].append(line)
            data['NonParametric_Flux'].append(res['flux_mean'][j])
            data['NonParametric_Flux_std'].append(res['flux_std'][j])
            data['Gaussian_Flux'].append(res['flux_fit_mean'][j])
            data['Gaussian_Flux_std'].append(res['flux_fit_std'][j])
            data['a'].append(res['a_mean'][j]); data['a_std'].append(res['a_std'][j])
            data['b'].append(res['b_mean'][j]); data['b_std'].append(res['b_std'][j])
            data['d'].append(res['d_mean'][j]); data['d_std'].append(res['d_std'][j])
            data['s'].append(res['s_mean'][j]); data['s_std'].append(res['s_std'][j])
            data['mu'].append(res['mu_mean'][j]); data['mu_std'].append(res['mu_std'][j])

            

    df = pd.DataFrame(data)
    df = df.replace(0.0, -9999).fillna(-9999) #filling in non-detections with -9999
    FILENAME = f'{output_dir}/combined_absorption_results.txt'
    filename = get_filename(output_dir, FILE_NAME)
    df.to_csv(filename, sep='\t', index=False)
    print(f"Results saved to: {filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit Gaussian absorption lines in observed frame with redshift table.")
    parser.add_argument("spectra_dir", help="Directory containing .txt spectra files")
    parser.add_argument("line_file", help="File containing line definitions and continuum windows")
    parser.add_argument("redshift_file", help="File containing object IDs and redshifts")
    parser.add_argument("output_dir", help="Directory to save results")

    args = parser.parse_args()
    main(args.spectra_dir, args.line_file, args.redshift_file, args.output_dir)
