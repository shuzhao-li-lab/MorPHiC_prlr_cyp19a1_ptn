<< comment
rm -rf lipidomics_differential_abundance_analysis
rm -rf annotated_significant_features
rm -rf mummichog_data
rm -rf figures
rm -rf mummichog_results

mkdir lipidomics_differential_abundance_analysis
mkdir annotated_significant_features
mkdir mummichog_data
mkdir figures
mkdir mummichog_results

python3 ./Analysis.py ./FeatureTables/Pellet_RP_Pos/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Pellet_Lipidomics+ /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_RP_Pos/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Pellets_18_Rppos Pellets_9_Rppos
python3 ./Analysis.py ./FeatureTables/Pellet_RP_Neg/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Pellet_Lipidomics- /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_RP_Neg/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Pellets_18_Rpneg Pellets_22_Rpneg
python3 ./Analysis.py ./FeatureTables/Pellet_HILIC_pos/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Pellet_Metabolomics+ /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_HILIC_pos/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Pellets_25_HILICpos Pellets_28_HILICpos
python3 ./Analysis.py ./FeatureTables/Pellet_HILIC_Neg/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Pellet_Metabolomics- /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_HILIC_Neg/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Pellets_25_HILICneg

python3 ./Analysis.py ./FeatureTables/Media_RP_Pos/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Supernatant_Lipidomics+ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_RP_Pos/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Media_7_Rppos Media_8_Rppos
python3 ./Analysis.py ./FeatureTables/Media_RP_Neg/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Supernatant_Lipidomics- /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_RP_Neg/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Media_18_Rpneg Media_7_Rpneg
python3 ./Analysis.py ./FeatureTables/Media_Hilic_Pos/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Supernatant_Metabolomics+ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_Hilic_Pos/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Media_18_HILICpos Media_29_HILICpos
python3 ./Analysis.py ./FeatureTables/Media_Hilic_Neg/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Supernatant_Metabolomics- /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_Hilic_Neg/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Media_6_HILICneg
comment
cd ./mummichog_results/
run_mummichog () {
    python3.8 -m mummichog.main -f /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/mummichog_data/$1 -o $2 -u 5
}

file_list=`/bin/ls /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/mummichog_data/ | grep -i metabolomics`
echo $file_list

for input_file in $file_list
do
    run_mummichog $input_file $input_file
done
cd ..