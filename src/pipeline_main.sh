#!/usr/bin/env bash

echo Running $(basename "${BASH_SOURCE}")

# Copy inputs to the working directory
copy_inputs.sh

# Convert eprime .txt to csv format
eprime_to_csv.py -o "${out_dir}"/eprime.csv "${out_dir}"/eprime.txt

# FSL based motion correction, topup, registration
fsl_processing.sh

# Unzip .nii ahead of SPM call
gunzip "${out_dir}"/ctrrfmri?.nii.gz \
    "${out_dir}"/ctrrfmri_mean_all.nii.gz \
    "${out_dir}"/biasnorm.nii.gz \
    "${out_dir}"/y_deffwd.nii.gz

# Handle missing fmri for SPM call
extra_cmd=
for r in 1 2 3 4; do
    if [[ ! -f "${out_dir}/ctrrfmri${r}.nii" ]]; then 
        extra_cmd="${extra_cmd} fmri${r}_nii NONE"
    fi
done

# SPM call
run_spm12.sh "${MATLAB_RUNTIME}" function matlab_entrypoint \
    hpf_sec "${hpf_sec}" \
    fwhm_mm "${fwhm_mm}" \
    out_dir "${out_dir}" \
    ${extra_cmd}

# Freeview-based PDF creation
make_pdf.sh

# Finalize and organize outputs
finalize.sh
