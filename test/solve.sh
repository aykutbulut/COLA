for input_file in r12c15k5i10.mps r14c18k3i18.mps r17c30k3i12.mps r22c40k10i30.mps  r23c45k3i27.mps  r27c50k5i25.mps r32c60k15i30.mps r52c75k5i40.mps
do
    echo "solving " $input_file
    #echo $input_file

    cola ../random/$input_file
    echo $input_file "done!"
done
