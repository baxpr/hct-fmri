#!/usr/bin/env bash
#
# Given a phase encoding direction as the first command line argument, generate
# a datain.txt file for topup. Assumes a two-volume nifti with forward and reverse
# phase encoding on the same axis, same effective echo spacing. No attempt made
# to use a correct value for the echo spacing.

pedir="${1}"

case "${pedir}" in

"+i")
cat > datain.txt <<HERE
1 0 0 1
-1 0 0 1
HERE
;;

"-i")
cat > datain.txt <<HERE
-1 0 0 1
1 0 0 1
HERE
;;

"+j")
cat > datain.txt <<HERE
0 1 0 1
0 -1 0 1
HERE
;;

"-j")
cat > datain.txt <<HERE
0 -1 0 1
0 1 0 1
HERE
;;

"+k")
cat > datain.txt <<HERE
0 0 1 1
0 0 -1 1
HERE
;;

"-k")
cat > datain.txt <<HERE
0 0 -1 1
0 0 1 1
HERE
;;

*)
echo Unknown phase encoding direction "${pedir}"
exit 1
;;

esac