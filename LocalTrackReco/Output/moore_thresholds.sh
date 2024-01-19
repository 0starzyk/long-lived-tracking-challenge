#!/bin/bash

algorithm_path="../MooreEnhanced/Pr/PrAlgorithms/src/PrLongLivedTracking.cpp"
baseline_path="../OptionFiles/Baseline.py"
moore_path="../MooreEnhanced/run"

for threshold in $(seq 0.01 0.01 0.99)
do
    sed -i -e "s/if (model<.*)continue;$/if (model<$threshold)continue;/" $algorithm_path
    sed -i -e "s/options\.histo_file = \"output_hist.*\.root\"/options\.histo_file = \"output_hist_$threshold\.root\"/g" $baseline_path
    sed -i -e "s/options\.ntuple_file = \"output_tuple.*\.root\"/options\.ntuple_file = \"output_tuple_$threshold\.root\"/g" $baseline_path
    make -C ../MooreEnhanced
    $moore_path gaudirun.py $baseline_path | tee logs/Enhance_$threshold.log
    mv Dumper_recTracks.root Dumper_recTracks_$threshold.root
done
