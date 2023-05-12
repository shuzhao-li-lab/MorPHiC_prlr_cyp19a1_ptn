<<comment
### FEATURE PROCESSING ###

# Pellet Lipidomics +
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Pellet_RP_Pos/ ./Sequence_iPSC_Pellets_regularmethods_03252023.csv pos --filter='{"Name": {"includes": ["rppos"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Pellet_RP_Pos/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Pellet_RP_Pos/experiment.json
#    The reference sample is:
#    ||* Pellets_45_Rppos *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Pos/experiment.json --table=preferred --all
# PCA - Samples cluster nicely apart from blanks and DDA. Pooled QC clusters with real samples.
# Correlation - looks good, samples correlate nicely with each other. Pooled correlate with all real samples
# Missing Feature Percentile - sigmoid
# Missing Feature Counts - very similar for real samples, previous outliers have more missing features
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Pellet_RP_Pos/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_RP_Pos/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Pos/experiment.json --table=processed_preferred --all --interactive
# Pellets 18 and Pellets 34 are outliers

# Pellet Metabolomics +
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Pellet_HILIC_pos/ ./Sequence_iPSC_Pellets_regularmethods_03252023.csv pos --filter='{"Name": {"includes": ["hilicpos"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Pellet_HILIC_pos/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Pellet_HILIC_pos/experiment.json
#     The reference sample is:
#    ||* Pellets_1_HILICpos *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=preferred --all
# PCA two clusters of samples, unsure which is which, complex vs basal media maybe?
# missing feature plot is not sigmoid, no clear inflection point.
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=processed_preferred --all --interactive
# Pellets_25_HILICpos, Pellets_5_HILICpos, and Pellets_28_HILICpos are outliers

# Pellet Lipidomics -
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Pellet_RP_Neg/ ./Sequence_iPSC_Pellets_regularmethods_03252023.csv neg --filter='{"Name": {"includes": ["rpneg"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Pellet_RP_Neg/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Pellet_RP_Neg/experiment.json
#    The reference sample is:
#    ||* Pellets_3_Rpneg *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Neg/experiment.json --table=preferred --all
# sigmoid, no clear PCA outliers in the samples, TSNE is a big blob
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Pellet_RP_Neg/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_RP_Neg/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Neg/experiment.json --table=processed_preferred --all --interactive
# Pellets_9_Rpneg, Pellets_39_Rpneg, Pellets_34_Rpneg are outliers in PCA, 39 and 34 have low correlation, 9 and 34 have many missing features too.

# Pellet Metabolomics - 
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Pellet_HILIC_Neg/ ./Sequence_iPSC_Pellets_regularmethods_03252023.csv neg --filter='{"Name": {"includes": ["hilicneg"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Pellet_HILIC_Neg/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Pellet_HILIC_Neg/experiment.json
#    The reference sample is:
#    ||* Pellets_9_HILICneg *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=preferred --all
# no pca outliers, sigmoid missing features
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=processed_preferred --all --interactive
# maybe pellets 11 is an outlier in TSNE?, Pellets_25_HILICneg and Pellets_30_HILICneg are outliers in missing features

python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/DmPA_pellets/ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_Pellets_TCA_Derivatization_03282023.csv pos
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/DmPA_pellets/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/DmPA_pellets/experiment.json
#    The reference sample is:
#    ||* Blank_TCA_derivatization_1 *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/DmPA_pellets/experiment.json --table=preferred --all --interactive
# sigmoid curve, decent looking correlation
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/DmPA_pellets/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/DmPA_pellets/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["ipsc"]}}' --drop_samples
# samples 13 and 14 failed to be converted it seems. Dropped manually from the table. 
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/DmPA_pellets/experiment.json --table=processed_preferred --all --interactive
# iPSC_Pellets_40_TCA_derivatization, iPSC_Pellets_54_TCA_derivatization, iPSC_Pellets_9_TCA_derivatization, and iPSC_Pellets_2_TCA_derivatization are outliers in PCA. 
# similar outliers in t-SNE.
# similar outliers in correlation heatmap
# these samples have many missing features
# 2, 54, and 9 are missing feature outliers


### ANNOTATION ###
# Lipidomics
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Pellet_RP_Pos/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Pellet_RP_Pos/experiment.json --table=default --new_table_moniker=LMSD_annotated_default ./LMSD.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Pellet_RP_Neg/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Pellet_RP_Neg/experiment.json --table=default --new_table_moniker=LMSD_annotated_default ./LMSD.json

# Metabolomics
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Pellet_HILIC_Neg/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=default --new_table_moniker=HMDB_annotated_default ./HMDB.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Pellet_HILIC_pos/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=default --new_table_moniker=HMDB_annotated_default ./HMDB.json
comment
### STATISTICAL ANALYSIS ###

python3 ./Analysis.py ./FeatureTables/Pellet_RP_Pos/filtered_feature_tables/processed_preferred_Feature_table.tsv  Pellet_Lipidomics+ /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_RP_Pos/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Pellets_18_Rppos Pellets_34_Rppos
python3 ./Analysis.py ./FeatureTables/Pellet_RP_Neg/filtered_feature_tables/processed_preferred_Feature_table.tsv  Pellet_Lipidomics- /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_RP_Neg/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Pellets_9_Rpneg Pellets_39_Rpneg Pellets_34_Rpneg
python3 ./Analysis.py ./FeatureTables/Pellet_HILIC_pos/filtered_feature_tables/processed_preferred_Feature_table.tsv  Pellet_Metabolomics+ /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_HILIC_pos/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Pellets_25_HILICpos Pellets_5_HILICpos Pellets_28_HILICpos
python3 ./Analysis.py ./FeatureTables/Pellet_HILIC_Neg/filtered_feature_tables/processed_preferred_Feature_table.tsv  Pellet_Metabolomics- /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_HILIC_Neg/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Pellets_25_HILICneg Pellets_30_HILICneg
