"""
=========================================================
bCFS Preprocessing in MNE-Python
kav bandara, unimelb 2025
=========================================================

This script performs Independent Component Analysis (ICA) on EEG/bCFS data. 

It loads the raw EEG data, applies preprocessing steps such as filtering and resampling,
and then runs ICA to identify and label components. The script generates visualizations of the ICA components
and saves the results as PDF reports and figures.

IClabel is used to automatically label the ICA components, and components that are likely to be artifacts
are excluded based on a specified threshold.
=========================================================
"""

import os.path as op
import os
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for HPC
import matplotlib.pyplot as plt
import shutil
import json
import argparse

from fpdf import FPDF
import mne
from mne.preprocessing import ICA
from mne_icalabel import label_components
import mne_bids

import sys
sys.path.insert(1, op.dirname(op.dirname(os.path.abspath(__file__))))

base_dir = r'/data/gpfs/projects/punim2118/bCFS/Data/EEG/'
# Get the path to the bCFS folder (parent of Data\EEG)
bcfs_root = os.path.dirname(os.path.dirname(base_dir))

# Define the path for the 'figures' and pdf folders
figure_root = os.path.join(bcfs_root, 'figures')
pdf_root = os.path.join(bcfs_root, 'reports')
os.makedirs(figure_root, exist_ok=True)
os.makedirs(pdf_root, exist_ok=True)

def load_raw(subject_id):
    
    """
    loads raw EEG data, resamples it to 200 hz, filters it at 1hz, and rereferences to the average

    parameters
    ----------
    subject_id : str
        The subject ID to load data for, e.g., "S09".

    returns
    -------
    raw : mne.io.Raw
        The loaded, resampled, and filtered MNE Raw object.
    """

    # Load the raw data
    print(f"--- Loading data for {subject_id} ---")
    fname = f"{subject_id}_bCFS.bdf" 
    fpath = op.join(base_dir, fname)

    # Load the file from the full path
    raw = mne.io.read_raw_bdf(fpath, preload=True)
    raw.load_data()

    # Correct channel labels 
    print("    -> Correcting channel types...")
    channel_type_mapping = {
        'LH': 'eog', 'RH': 'eog',
        'LV': 'eog', 'UV': 'eog',
        'M1': 'misc', 'M2': 'misc',
        'Nz': 'misc', 'SNz': 'misc'
    }
    raw.set_channel_types(mapping=channel_type_mapping)
    
    # set template montage
    try:
        raw.set_montage('standard_1020')
    except ValueError as e:

    raw.resample(200, npad="auto") 
    #raw.filter(l_freq=1.0, h_freq=None)
    #raw.set_eeg_reference(ref_channels='average')  # Set average reference

    return raw 

def run_ica(subject_id, threshold=0.95):

    """

    Generates ICA decompositions for EEG, plots
    component timecourses and topographies, saves ICA solutions, and
    generates PDF reports. 
    
    Includes automatic ICA component labeling using MNE-ICLabel.

    Parameters
    ----------
    subject_id : str
        Subject ID.

    threshold : float, optional
        Threshold for automatic component exclusion. Defaults to 0.95.

    """

    raw = load_raw(subject_id)
    
    raw_for_ica_fit = raw.copy().filter(l_freq=1.0, h_freq=None)
    
    ica = ICA(
        n_components=0.99, 
        max_iter="auto",
        method="infomax",
        random_state=97,
        fit_params=dict(extended=True))

    # Define rejection criteria to exclude noisy segments during fitting -- this is because if the data is excessively noisy the ica algorithm cant find many components to explain the data
    reject_criteria = dict(eeg=250e-6)  # 200 µV
    
    # Run ICA on filtered raw data    
    ica.fit(raw_for_ica_fit,
            reject = reject_criteria,
            picks='eeg',
            verbose=True)
    
    # =========================
    # Selecting ICA components automatically
    # =========================
    
    # Applies the automatic ICA component labeling algorithm, which will
    # assign a probability value for each component being one of:
    
    # - brain
    # - muscle artifact
    # - eye blink
    # - heart beat
    # - line noise
    # - channel noise
    # - other
 
    # prepare a temporary, average-referenced copy of the data for ICLabel.
    raw_for_iclabel = raw_for_ica_fit.copy().set_eeg_reference('average', projection=False)
    ic_labels = label_components(raw_for_iclabel, ica, method="iclabel")
    labels = ic_labels["labels"]
    prob = ic_labels["y_pred_proba"]
    ica.labels_ = labels # Store labels in ICA object for report later

    # find and exclude components based on the labels and probabilities
    ica_exclude_idx = []
    threshold = threshold  # Set a threshold for component exclusion

    for i, label in enumerate(ic_labels['labels']):
        # We only care about specific artifact labels when the model is 95% confident
        # that the component is an artifact.
        if label in ['muscle artifact', 'eye blink', 'heart beat', 'line noise', 'channel noise']:
            probability = ic_labels['y_pred_proba'][i]
            if probability > threshold:
                ica_exclude_idx.append(i)

    print(f"Using a threshold of {threshold}, suggested components to exclude: {ica_exclude_idx}")

    # Save ICA file
    fname_ica = op.join(base_dir, f'{subject_id}-ica.fif')
    ica.save(fname_ica, overwrite=True) # Save ICA object to the new path
    
    # =========================
    # Save Exclusion Information
    # =========================

    # Save to JSON file
    excludes_json = op.join(base_dir, "exclude_ics.json") 

    try:
        with open(excludes_json, 'r') as f:
            all_ica_rej = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        # If the file doesn't exist yet, or is empty, start a new dictionary
        print(f"Master JSON not found or empty. Creating a new one at: {excludes_json}")
        all_ica_rej = {}

    this_sub_ica_rej = {
        subject_id: {
                "eeg_ica_excluded": ica_exclude_idx,
        }
    }
        
    all_ica_rej[subject_id] = this_sub_ica_rej

    with open(excludes_json, 'w') as f:
        json.dump(all_ica_rej, f, indent=4)

    #Generate pdf for this subject
    pdf = FPDF(orientation="P", unit="mm", format="A4") 
    effective_width = pdf.w - pdf.l_margin - pdf.r_margin

    # helper function to handle single figure or list of figures
    def save_and_add_figs_to_pdf(fig_object, base_fname, title_prefix, pdf_obj):
        if not isinstance(fig_object, list):
            fig_object = [fig_object]  # If it's a single fig, put it in a list
        
        for i, fig in enumerate(fig_object):
            fname = op.join(figure_root, f"{base_fname}_part{i}.png")
            fig.savefig(fname, dpi=150)
            plt.close(fig)
            
            pdf_obj.add_page()
            pdf_obj.set_font('helvetica', 'B', 16)
            pdf_obj.cell(0, 10, f"{title_prefix} (Part {i+1})", ln=1)
            pdf_obj.image(fname, x=pdf_obj.l_margin, y=25, w=effective_width)

    # Plot timecourse and topography of the ICs in n sets of 20
    n_comp_list = range(ica.n_components_)
    plot_comp_list = [n_comp_list[i:i + 20] for i in range(0, len(n_comp_list), 20)]
    
    # Plot overview
    figs_overview = ica.plot_components(show=False)
    save_and_add_figs_to_pdf(figs_overview, f"{subject_id}_ica_overview", f"{subject_id} - IC Topographies", pdf)

    # Plot time courses
    n_comp_list = range(ica.n_components_)
    plot_comp_list = [n_comp_list[i:i + 20] for i in range(0, len(n_comp_list), 20)]
    for i, comps in enumerate(plot_comp_list):
        fig_src = ica.plot_sources(raw, picks=comps, show_scrollbars=False, show=False)
        save_and_add_figs_to_pdf(fig_src, f"{subject_id}_ica_sources_batch{i}", f"{subject_id} - IC Time Courses (ICs {comps[0]}-{comps[-1]})", pdf)

    # Add probability table
    pdf.add_page()
    pdf.set_font('helvetica', 'B', 16)
    pdf.cell(0, 10, f"{subject_id} - ICLabel Probabilities", ln=1)
    pdf.ln(5)
    pdf.set_font('Courier', '', 10) 
    table_header = f"{'IC':<5} | {'Predicted Label':<20} | {'Confidence':<10}\n"
    table_header += "-" * len(table_header) + "\n"
    table_rows = ""
    for idx, label in enumerate(ic_labels['labels']):
        prob = ic_labels['y_pred_proba'][idx]
        is_excluded = "*" if idx in ica.exclude else ""
        table_rows += f"{idx:<5} | {label:<20} | {prob:<10.4f}{is_excluded}\n"
    pdf.multi_cell(0, 5, table_header + table_rows)
    pdf.set_font('helvetica', 'I', 8)
    pdf.cell(0, 5, "* denotes components marked for exclusion.", ln=1)
    
    fname_report = op.join(pdf_root, f"{subject_id}_ICA-report.pdf")
    pdf.output(fname_report)
    print(f"--- Report saved to {fname_report} ---")  

    # applies ICA solution

    print(f"--- Applying ICA solution for {subject_id} ---")
    print(f"--- Removing ICA components: {ica_exclude_idx} ---")

    ica.apply(raw, exclude=ica_exclude_idx)
    
    # apply the average reference to the cleaned data
    print("--- Applying average EEG reference to cleaned data ---")
    raw.set_eeg_reference(ref_channels='average', projection = False)

    # Save filtered data
    fname_filt = f"{subject_id}_bCFS_ica_filt.fif"
    fpath_filt = op.join(base_dir, fname_filt)
    
    raw.save(fpath_filt, overwrite=True)
    
#execute script on HPC
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run ICA pipeline for a single subject.")
    parser.add_argument('--subject', type=str, required=True, help='The subject ID to process (e.g., S26).')
    args = parser.parse_args()
    
    print(f"--- Starting full ICA pipeline for subject: {args.subject} ---")
    run_ica(args.subject, threshold=0.95)
    print(f"\n--- SCRIPT COMPLETE for subject {args.subject} ---")
                 
"""
if __name__ == '__main__':
    
    # SPECIFY THE SUBJECT TO PROCESS HERE
    # ===================================================================
    subject_to_process = "S26" 
    # ===================================================================
    
    print(f"--- Starting ica preprocessing for subject: {subject_to_process} ---")
    
    run_ica(subject_to_process)

    print(f"--- Finished ica preprocessing for subject: {subject_to_process} ---")
    
"""
