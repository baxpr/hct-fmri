#!/usr/bin/env bash
#
# Take two single fmri volumes as input and run topup. Filenames must be
# given basename only and the "FSL way" without extension. Runs in the current
# directory. See make_datain.sh for pedir options.
#
# run_topup.sh <pedir> <fwd_meanfmri> <rev_meanfmri> <fwd_timeseries_fmri>

pedir="${1}"
fwd="${2}"
rev="${3}"
run1="${4}"
run2="${5}"
run3="${6}"
run4="${7}"

# Combine fMRIs to make topup input with two vols
fslmerge -t topup "${fwd}" "${rev}"

# Make topup phase dir info file
make_datain.sh "${pedir}"

# Run topup. Save the field map, but realize the sign and amplitude are not meaningful.
# Use b02b0_1 schedule to avoid issue with odd number of slices
topup --imain=topup --datain=datain.txt --fout=topup_fieldmap_hz --config=b02b0_1.cnf

# Apply topup to fwd (this is the data we will actually analyze for fmri)
applytopup --imain="${fwd}" --inindex=1 --datain=datain.txt \
	--topup=topup --out=t"${fwd}" --method=jac

# Apply topup to rev - just for reference
applytopup --imain="${rev}" --inindex=2 --datain=datain.txt \
	--topup=topup --out=t"${rev}" --method=jac

# Apply topup to actual time series
for f in "${run1}" "${run2}" "${run3}" "${run4}"; do
    applytopup --imain="${f}" --inindex=1 --datain=datain.txt \
	    --topup=topup --out=t"${f}" --method=jac
done
