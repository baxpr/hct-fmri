#!/usr/bin/env bash

echo Running $(basename "${BASH_SOURCE}")

cd "${out_dir}"

# Zip nifti files in SPM outputs
for d in \
    spm_cue_orth0 \
    spm_cue_orth1_epsi2 \
    spm_cue_orth1_mu33 \
    spm_cue_mu3_feedback_epsi3_orth1 \
; do
    gzip "${d}"/*.nii
done

# Zip unsmoothed mean fmri
gzip wctrrfmri_mean_all.nii

# Preprocessed fmri
mkdir SWFMRI
cp swctrrfmri?.nii SWFMRI
gzip SWFMRI/*.nii

