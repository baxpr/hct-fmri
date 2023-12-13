#!/usr/bin/env bash

echo Making PDF

# Work in output directory
cd ${out_dir}

# Timestamp
thedate=$(date)


# EPI-to-T1 coregistration, native space, near WM COM
com=( $(fslstats wm -c) )
com[0]=$(echo ${com[0]} + 10 | bc)
fsleyes render -of coreg.png \
	--scene ortho --worldLoc ${com[@]} --displaySpace world --size 1800 600 --xzoom 1000 --yzoom 1000 --zzoom 1000 \
	--layout horizontal --hideCursor \
	ctrrfmri_mean_all --overlayType volume \
	gm --overlayType label --outline --outlineWidth 2 --lut harvard-oxford-subcortical


# EPI normalization, MNI space
fslmaths ${FSLDIR}/data/standard/tissuepriors/avg152T1_gray -thr 100 -bin gm_mni
fsleyes render -of mni.png \
	--scene ortho --worldLoc 10 -20 0 --displaySpace world --size 1800 600 --xzoom 600 --yzoom 600 --zzoom 600 \
	--layout horizontal --hideCursor \
	wctrrfmri_mean_all --overlayType volume \
	gm_mni --overlayType label --outline --outlineWidth 2 --lut harvard-oxford-subcortical


# fMRI contrast image, slices
spm_dir=spm_cue_orth1_epsi2
for connum in 1; do
    connum0=$(printf "%04g\n" ${connum})
    conname=$(get_conname.py ${out_dir}/spm_contrast_names_cue_orth1_epsi2.csv ${connum})
    c=10
    for slice in -35 -20 -5 10 25 40 55 70  ; do
	    ((c++))
	    fsleyes render -of ${spm_dir}_${connum0}_${c}.png \
	        --scene ortho --worldLoc 0 0 ${slice} --displaySpace world --size 600 600 --yzoom 1000 \
	        --layout horizontal --hideCursor --hideLabels --hidex --hidey \
		    biasnorm --overlayType volume \
		    ${spm_dir}/spmT_${connum0} --overlayType volume --displayRange 3 10 \
		    --useNegativeCmap --cmap red-yellow --negativeCmap blue-lightblue
    done

    montage \
	    -mode concatenate ${spm_dir}_${connum0}_??.png \
	    -tile 3x -quality 100 -background black -gravity center \
	    -border 20 -bordercolor black ${spm_dir}_${connum0}.png

    convert -size 2600x3365 xc:white \
	    -gravity center \( ${spm_dir}_${connum0}.png -resize 2400x \) -composite \
	    -gravity North -pointsize 48 -annotate +0+100 \
	    "PRL fMRI, ${spm_dir}, contrast ${connum}: ${conname}" \
	    -gravity SouthEast -pointsize 48 -annotate +100+100 "${thedate}" \
	    page_${spm_dir}_${connum0}.png

done


# Combine
convert -size 2600x3365 xc:white \
	-gravity center \( coreg.png -geometry '2400x2400+0-0' \) -composite \
	-gravity center -pointsize 48 -annotate +0-1250 \
	"T1 gray matter outline on unsmoothed registered mean fMRI (native space)" \
	-gravity center \( mni.png -geometry '2400x2400+0+1200' \) -composite \
	-gravity center -pointsize 48 -annotate +0-50 \
	"Atlas gray matter outline on unsmoothed warped mean fMRI (atlas space)" \
	-gravity SouthEast -pointsize 48 -annotate +100+100 "${thedate}" \
	page_reg.png

convert -size 2600x3365 xc:white \
	-gravity center \( first_level_design_cue_orth1_epsi2_001.png -resize 2000x \) -composite \
	-gravity SouthEast -pointsize 48 -annotate +100+100 "${thedate}" \
	page_design_cue_orth1_epsi2.png

convert \
    page_reg.png \
    page_design_cue_orth1_epsi2.png page_spm_cue_orth1_epsi2_*.png \
    prl-fmri.pdf

