#!/usr/bin/env bash
#
# Motion correction, topup, and registration to T1 for 4 fMRI time series 
# with a matched time series or volume acquired with reverse phase encoding
# direction.
#
# Relies on env vars exported from pipeline_entrypoint.sh to get arguments:
#    out_dir
#    pedir
#    vox_mm
#
# Assumed filenames in $out_dir are
#    seg.nii.gz          Multiatlas/slant segmentation result
#    biascorr.nii.gz     CAT12 bias-corrected T1 in native space
#    icv.nii.gz          T1 masked to brain only
#    fmri?.nii.gz        fMRI time series 1-4
#    fmritopup.nii.gz    fMRI with reversed phase encoding direction
#
# Results are
#    ctrrfmri_mean_all.nii.gz        Mean of fMRIs, topup'd and registered to T1
#    ctrrfmri?.nii.gz                Topup'd and registered fMRI time series
#    ctrrfmritopup_mean_reg.nii.gz   Mean rev phase enc fMRI after topup/reg

# Get in working dir
cd "${out_dir}"

# Gray matter mask from slant
fslmaths seg -thr  3.5 -uthr  4.5 -bin -mul -1 -add 1 -mul seg tmp
fslmaths tmp -thr 10.5 -uthr 11.5 -bin -mul -1 -add 1 -mul tmp tmp
fslmaths tmp -thr 39.5 -uthr 41.5 -bin -mul -1 -add 1 -mul tmp tmp
fslmaths tmp -thr 43.5 -uthr 45.5 -bin -mul -1 -add 1 -mul tmp tmp
fslmaths tmp -thr 48.5 -uthr 52.5 -bin -mul -1 -add 1 -mul tmp tmp
fslmaths tmp -bin gm
rm tmp.nii.gz

# White matter mask from slant
echo "White matter mask"
fslmaths seg -thr 39.5 -uthr 41.5 -bin tmp
fslmaths seg -thr 43.5 -uthr 45.5 -add tmp -bin wm
rm tmp.nii.gz

# Motion correction within run, and for the short topup series
echo "Motion correction"
for n in 1 2 3 4; do
    echo "    Run ${n}"
    mcflirt -in fmri${n} -meanvol -out rfmri${n} -plots
done

echo "    Topup run"
mcflirt -in fmritopup -meanvol -out rfmritopup


# Alignment between runs and overall mean fmri
echo "Aligning runs"
cp rfmri1_mean_reg.nii.gz rrfmri1_mean_reg.nii.gz
cp rfmri1.nii.gz rrfmri1.nii.gz
opts="-usesqform -searchrx -5 5 -searchry -5 5 -searchrz -5 5"
for n in 2 3 4; do
    echo "    Run ${n} to run 1"
    flirt ${opts} -in rfmri${n}_mean_reg -ref rrfmri1_mean_reg \
        -out rrfmri${n}_mean_reg -omat r${n}to1.fslmat
    flirt -applyxfm -init r${n}to1.fslmat -in rfmri${n} -ref rrfmri1_mean_reg -out rrfmri${n}
done

echo "    Topup run to run 1"
flirt ${opts} -in rfmritopup_mean_reg -ref rrfmri1_mean_reg -out rrfmritopup_mean_reg

echo "    Computing overall mean"
fslmaths rrfmri1_mean_reg -add rrfmri2_mean_reg \
    -add rrfmri3_mean_reg -add rrfmri4_mean_reg \
    -div 4 rrfmri_mean_all


# Run topup. After this, the 'tr' prefix files always contain the data that will be further
# processed.
echo "Running TOPUP"
run_topup.sh "${pedir}" rrfmri_mean_all rrfmritopup_mean_reg rrfmri1 rrfmri2 rrfmri3 rrfmri4

# Register corrected mean fmri to T1. biascorr is the adjusted T1 from cat12, ICV is the 
# ICV_NATIVE resource from cat12 that is masked to only brain.
echo "Coregistration"
epi_reg --epi=trrfmri_mean_all --t1=biascorr --t1brain=icv --wmseg=wm --out=ctrrfmri_mean_all
mv ctrrfmri_mean_all.mat ctrrfmri_mean_all.fslmat

# Use flirt to resample to the desired voxel size, overwriting epi_reg output image
flirt -applyisoxfm "${vox_mm}" -init ctrrfmri_mean_all.fslmat -in trrfmri_mean_all \
	-ref biascorr -out ctrrfmri_mean_all

# Apply coregistration to the corrected time series
for n in 1 2 3 4; do
    flirt -applyisoxfm "${vox_mm}" -init ctrrfmri_mean_all.fslmat \
        -in trrfmri${n} -ref biascorr -out ctrrfmri${n}
done

# And to the topup image, for reference
flirt -applyisoxfm "${vox_mm}" -init ctrrfmri_mean_all.fslmat \
    -in trrfmritopup_mean_reg -ref biascorr -out ctrrfmritopup_mean_reg

