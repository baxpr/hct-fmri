#!/usr/bin/env bash

echo Running $(basename "${BASH_SOURCE}")

# Files
if [[ -f "${fmri1_niigz}" ]]; then
    cp "${fmri1_niigz}" "${out_dir}"/fmri1.nii.gz
fi
if [[ -f "${fmri2_niigz}" ]]; then
    cp "${fmri2_niigz}" "${out_dir}"/fmri2.nii.gz
fi
if [[ -f "${fmri3_niigz}" ]]; then
    cp "${fmri3_niigz}" "${out_dir}"/fmri3.nii.gz
fi
if [[ -f "${fmri4_niigz}" ]]; then
    cp "${fmri4_niigz}" "${out_dir}"/fmri4.nii.gz
fi
cp "${fmritopup_niigz}" "${out_dir}"/fmritopup.nii.gz
cp "${seg_niigz}" "${out_dir}"/seg.nii.gz
cp "${icv_niigz}" "${out_dir}"/icv.nii.gz
cp "${deffwd_niigz}" "${out_dir}"/y_deffwd.nii.gz
cp "${biascorr_niigz}" "${out_dir}"/biascorr.nii.gz
cp "${biasnorm_niigz}" "${out_dir}"/biasnorm.nii.gz
cp "${eprime_txt}" "${out_dir}"/eprime.txt
