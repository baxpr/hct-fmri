#!/usr/bin/env bash

echo Running $(basename "${BASH_SOURCE}")

cd "${out_dir}"

# Zip nifti files in SPM outputs
for d in \
    spm_hct \
; do
    gzip "${d}"/*.nii
done

# Zip unsmoothed mean fmri
gzip wctrrfmri_mean_all.nii

# Preprocessed fmri
mkdir SWFMRI
cp swctrrfmri?.nii SWFMRI
gzip SWFMRI/*.nii

