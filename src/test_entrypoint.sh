#!/usr/bin/env bash

export PATH=$(pwd):$(pwd)/../matlab/bin:$PATH

pipeline_entrypoint.sh \
    --fmri1_niigz $(pwd)/../INPUTS/run1.nii.gz \
    --fmri2_niigz $(pwd)/../INPUTS/run2.nii.gz \
    --fmri3_niigz $(pwd)/../INPUTS/run3.nii.gz \
    --fmri4_niigz NONE \
    --fmritopup_niigz $(pwd)/../INPUTS/topup.nii.gz \
    --seg_niigz $(pwd)/../INPUTS/seg.nii.gz \
    --icv_niigz $(pwd)/../INPUTS/p0t1.nii.gz \
    --refimg_nii avg152T1.nii \
    --deffwd_niigz $(pwd)/../INPUTS/y_t1.nii.gz \
    --biascorr_niigz $(pwd)/../INPUTS/mt1.nii.gz \
    --biasnorm_niigz $(pwd)/../INPUTS/wmt1.nii.gz \
    --eprime_txt $(pwd)/../INPUTS/eprime.txt \
    --pedir '+j' \
    --vox_mm 2 \
    --hpf_sec 200 \
    --fwhm_mm 6 \
    --out_dir $(pwd)/../OUTPUTS

