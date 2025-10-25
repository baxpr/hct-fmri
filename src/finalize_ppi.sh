#!/usr/bin/env bash

echo Running $(basename "${BASH_SOURCE}")

cd "${out_dir}"

# Clean up spm_hct
rm spm_hct/PPI*
rm spm_hct/VOI*

# Zip nifti files in SPM outputs
for d in spm_hct spm_hct_manppi* ppi; do
    gzip "${d}"/*.nii
done

# Zip unsmoothed mean fmri
gzip wctrrfmri_mean_all.nii

# Preprocessed fmri
mkdir WFMRI
cp wctrrfmri?.nii WFMRI
gzip WFMRI/*.nii
mkdir SWFMRI
cp swctrrfmri?.nii SWFMRI
gzip SWFMRI/*.nii

# In yaml, capture one ROI result in full as an example: spm_hct_manppi_dAI_L

# Rename PPI/VOI dir
mv ppi PPI_VOI_WORKFILES

# Grab just the PPI con images for small and coherent PPI output. Note exact
# specific assumption about which con image is which in the first level stats code
mkdir PPI_RESULTS
for roi in dAI_L dAI_R vAI_L vAI_R PI_L PI_R dAmy_L dAmy_R aMCC_L aMCC_R pACC_L pACC_R sgACC_L sgACC_R; do
    cp spm_hct_manppi_${roi}/con_0005.nii.gz \
        PPI_RESULTS/ppicon_0005_HeartConn_${roi}.nii.gz
    cp spm_hct_manppi_${roi}/con_0006.nii.gz \
        PPI_RESULTS/ppicon_0006_CountingConn_${roi}.nii.gz
    cp spm_hct_manppi_${roi}/con_0007.nii.gz \
        PPI_RESULTS/ppicon_0007_HeartConn_plus_CountingConn_${roi}.nii.gz
    cp spm_hct_manppi_${roi}/con_0008.nii.gz \
        PPI_RESULTS/ppicon_0008_HeartConn_gt_CountingConn_${roi}.nii.gz
done

# Capture the atlas
mkdir ATLAS
cp $(dirname "${BASH_SOURCE}")/matlab/src/atlas* ATLAS

