rm -rf lipidomics_differential_abundance_analysis
rm -rf annotated_significant_features
rm -rf mummichog_data
rm -rf significant_features
rm -rf FeatureTables
rm -rf figures
mkdir lipidomics_differential_abundance_analysis
mkdir annotated_significant_features
mkdir mummichog_data
mkdir FeatureTables
mkdir figures
sleep 1

### PELLET PROCESSING ###

# Pellet Lipidomics +
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Pellet_RP_Pos/ ./Sequence_iPSC_Pellets_regularmethods_03252023.csv pos --filter='{"Name": {"includes": ["rppos"], "lacks": ["qstd", "dda"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Pellet_RP_Pos/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
sleep 5
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
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Pos/experiment.json --table=processed_preferred --all
# Pellets 9 are outlier in PCA and missing feature
# Pellets 18 outlier in missing features
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_RP_Pos/experiment.json --table=preferred --new_table_moniker=processed_gfiltering .80 0 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Pos/experiment.json --table=processed_gfiltering --all
# no outliers

# Pellet Metabolomics +
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Pellet_HILIC_pos/ ./Sequence_iPSC_Pellets_regularmethods_03252023.csv pos  --filter='{"Name": {"includes": ["hilicpos"], "lacks": ["qstd", "dda"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Pellet_HILIC_pos/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
sleep 5
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Pellet_HILIC_pos/experiment.json
#     The reference sample is:
#    ||* Pellets_1_HILICpos *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=preferred --all
# PCA two clusters of samples, unsure which is which, complex vs basal media maybe?
# missing feature plot is not sigmoid, no clear inflection point.
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=processed_preferred --all
# Pellets_25_HILICpos and Pellets_28_HILICpos are outliers in PCA and missing features
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=preferred --new_table_moniker=processed_gfiltering .80 0 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=processed_gfiltering --all
# No outliers

# Pellet Lipidomics -
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Pellet_RP_Neg/ ./Sequence_iPSC_Pellets_regularmethods_03252023.csv neg  --filter='{"Name": {"includes": ["rpneg"], "lacks": ["qstd", "dda"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Pellet_RP_Neg/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
sleep 5
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Pellet_RP_Neg/experiment.json
#    The reference sample is:
#    ||* Pellets_3_Rpneg *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Neg/experiment.json --table=preferred --all
# sigmoid, no clear PCA outliers in the samples, TSNE is a big blob
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Pellet_RP_Neg/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_RP_Neg/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Neg/experiment.json --table=processed_preferred --all
# Pellets_22_Rpneg, Pellets_39_Rpneg, Pellets_54_Rpneg are outliers in PCA, 39 and 34 have low correlation, 9 and 34 have many missing features too.
# Pellets 18 and 22 are outliers in missing features
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_RP_Neg/experiment.json --table=preferred --new_table_moniker=processed_gfiltering .80 0 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_RP_Neg/experiment.json --table=processed_gfiltering --all
# Pellets 3 is an outlier in the sense that it has a Z of > -2.5 in missing feature count but is that a bad thing?

# Pellet Metabolomics - 
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Pellet_HILIC_Neg/ ./Sequence_iPSC_Pellets_regularmethods_03252023.csv neg  --filter='{"Name": {"includes": ["hilicneg"], "lacks": ["qstd", "dda"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Pellet_HILIC_Neg/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
sleep 5
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Pellet_HILIC_Neg/experiment.json
#    The reference sample is:
#    ||* Pellets_9_HILICneg *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=preferred --all
# no pca outliers, sigmoid missing features
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=processed_preferred --all
# maybe pellets 54 is an outlier in TSNE?
# Pellets_25_HILICneg is an outlier in missing feature counts
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=preferred --new_table_moniker=processed_gfiltering .80 0 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["pellets"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=processed_gfiltering --all


### SUPERNATANT PROCESSING###

# Supernatant lipidomics + 
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Media_RP_Pos/ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv pos --filter='{"Name": {"includes": ["rppos"], "lacks": ["qstd", "dda"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Media_RP_Pos/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
sleep 5
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Media_RP_Pos/experiment.json
#    The reference sample is:
#    ||* Pool_QC_Media_RPpos_3 *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_RP_Pos/experiment.json --table=preferred --all
# Media 12 might be an outlier in PCA, but not a crazy outlier. Relatively low correlation as well. 
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Media_RP_Pos/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Media_RP_Pos/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["media"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_RP_Pos/experiment.json --table=processed_preferred --all
# Media 6 and 47 are outliers in TSNE, but no obvious outlier in PCA
# Media 7 and 8 are missing feature outliers
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Media_RP_Pos/experiment.json --table=preferred --new_table_moniker=processed_gfiltering .80 0 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["media"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_RP_Pos/experiment.json --table=processed_gfiltering --all
# Media 12 is major outlier in PCA and is *almost* an outlier in median correlation z-score
# Media 28 is an outlier in missing feature but with more features

# Supernantant lipidomics -
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Media_RP_Neg/ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv neg  --filter='{"Name": {"includes": ["rpneg"], "lacks": ["qstd", "dda"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Media_RP_Neg/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
sleep 5
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Media_RP_Neg/experiment.json
#    The reference sample is:
#    ||* Pool_QC_Media_RPneg_3 *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_RP_Neg/experiment.json --table=preferred --all
# media 12 might be outlier in PCA, media 21 as well. Media 12 is an outlier in TSNE. 12 and 21 low correlation (relative)
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Media_RP_Neg/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Media_RP_Neg/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["media"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_RP_Neg/experiment.json --table=processed_preferred --all
# Media 18 and Media 7 are outliers in missing features. 
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Media_RP_Neg/experiment.json --table=preferred --new_table_moniker=processed_gfiltering .80 0 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["media"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_RP_Neg/experiment.json --table=processed_gfiltering --all
# Media 12 is outlier in PCA but this PCA looks weird, many clusters. 
# Media 21 is an outlier in missing feature but with more features

# Supernatant Metabolomics +
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Media_Hilic_Pos/ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv pos  --filter='{"Name": {"includes": ["hilicpos"], "lacks": ["qstd", "dda"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Media_Hilic_Pos/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
sleep 5
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Media_Hilic_Pos/experiment.json
#    The reference sample is:
#    ||* Media_14_HILICpos *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_Hilic_Pos/experiment.json --table=preferred --all 
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Media_Hilic_Pos/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Media_Hilic_Pos/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["media"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_Hilic_Pos/experiment.json --table=processed_preferred --all 
# Media_45_HILICpos PCA / TSNE outliers
# Media_18_HILICpos, Media_29_HILICpos are outliers in missing feature
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Media_Hilic_Pos/experiment.json --table=preferred --new_table_moniker=processed_gfiltering .80 0 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["media"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_Hilic_Pos/experiment.json --table=processed_gfiltering --all 
# Media 42 maybe PCA outlier, but again, multiple clusters. 

# Supernatant Metabolomics - 
python3 ../PythonCentricPipelineForMetabolomics/src/main.py assemble_experiment_from_CSV ./FeatureTables/Media_Hilic_Neg/ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv neg  --filter='{"Name": {"includes": ["hilicneg"], "lacks": ["qstd", "dda"]}}'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py convert_to_mzML ./FeatureTables/Media_Hilic_Neg/experiment.json $(which mono) ../PythonCentricPipelineForMetabolomics/ThermoRawFileParser/ThermoRawFileParser.exe
sleep 5
python3 ../PythonCentricPipelineForMetabolomics/src/main.py asari_full_processing ./FeatureTables/Media_Hilic_Neg/experiment.json
#    The reference sample is:
#    ||* Media_50_HILICneg *||
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_Hilic_Neg/experiment.json --table=preferred --all
python3 ../PythonCentricPipelineForMetabolomics/src/main.py drop_samples ./FeatureTables/Media_Hilic_Neg/experiment.json --table=preferred --substring_name_drop='["blank", "pool", "qstd"]'
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Media_Hilic_Neg/experiment.json --table=preferred --new_table_moniker=processed_preferred .80 .80 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["media"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_Hilic_Neg/experiment.json --table=processed_preferred --all 
# Meia_6_HILICneg is an outlier in missing features, maybe an outlier in PCA
python3 ../PythonCentricPipelineForMetabolomics/src/main.py preprocess_features ./FeatureTables/Media_Hilic_Neg/experiment.json --table=preferred --new_table_moniker=processed_gfiltering .80 0 3 '{"Name": {"includes": ["blank"]}}' '{"Name": {"includes": ["media"]}}' --drop_samples
python3 ../PythonCentricPipelineForMetabolomics/src/main.py feature_QCQA ./FeatureTables/Media_Hilic_Neg/experiment.json --table=processed_gfiltering --all
# Media 42 is an outlier maybe but also multiple clusters in PCA
# Media 6 is an outlier in missing features!

### ANNOTATION ###
# Pellet Lipidomics
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Pellet_RP_Pos/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Pellet_RP_Pos/experiment.json --table=default --new_table_moniker=LMSD_annotated_default ./LMSD.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Pellet_RP_Neg/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Pellet_RP_Neg/experiment.json --table=default --new_table_moniker=LMSD_annotated_default ./LMSD.json

# Pellet Metabolomics
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Pellet_HILIC_Neg/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Pellet_HILIC_Neg/experiment.json --table=default --new_table_moniker=HMDB_annotated_default ./HMDB.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Pellet_HILIC_pos/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Pellet_HILIC_pos/experiment.json --table=default --new_table_moniker=HMDB_annotated_default ./HMDB.json


### ANNOTATION ###
# Supernatant Lipidomics
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Media_RP_Pos/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Media_RP_Pos/experiment.json --table=default --new_table_moniker=LMSD_annotated_default ./LMSD.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Media_RP_Neg/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Media_RP_Neg/experiment.json --table=default --new_table_moniker=LMSD_annotated_default ./LMSD.json

# Supernatant Metabolomics
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Media_Hilic_Pos/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Media_Hilic_Pos/experiment.json --table=default --new_table_moniker=HMDB_annotated_default ./HMDB.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py build_empCpds ./FeatureTables/Media_Hilic_Neg/experiment.json
python3 ../PythonCentricPipelineForMetabolomics/src/main.py MS1_annotate ./FeatureTables/Media_Hilic_Neg/experiment.json --table=default --new_table_moniker=HMDB_annotated_default ./HMDB.json


### STATISTICAL ANALYSIS ###
# Most filtered tables #
#python3 ./Analysis.py ./FeatureTables/Pellet_RP_Pos/filtered_feature_tables/processed_preferred_Feature_table.tsv  Pellet_Lipidomics+ /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_RP_Pos/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Pellets_18_Rppos Pellets_9_Rppos
#python3 ./Analysis.py ./FeatureTables/Pellet_RP_Neg/filtered_feature_tables/processed_preferred_Feature_table.tsv  Pellet_Lipidomics- /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_RP_Neg/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Pellets_18_Rpneg Pellets_22_Rpneg
#python3 ./Analysis.py ./FeatureTables/Pellet_HILIC_pos/filtered_feature_tables/processed_preferred_Feature_table.tsv  Pellet_Metabolomics+ /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_HILIC_pos/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Pellets_25_HILICpos Pellets_28_HILICpos
#python3 ./Analysis.py ./FeatureTables/Pellet_HILIC_Neg/filtered_feature_tables/processed_preferred_Feature_table.tsv  Pellet_Metabolomics- /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_HILIC_Neg/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Pellets_25_HILICneg

#python3 ./Analysis.py ./FeatureTables/Media_RP_Pos/filtered_feature_tables/processed_preferred_Feature_table.tsv  Supernatant_Lipidomics+ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_RP_Pos/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Media_7_Rppos Media_8_Rppos
#python3 ./Analysis.py ./FeatureTables/Media_RP_Neg/filtered_feature_tables/processed_preferred_Feature_table.tsv  Supernatant_Lipidomics- /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_RP_Neg/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Media_18_Rpneg Media_7_Rpneg
#python3 ./Analysis.py ./FeatureTables/Media_Hilic_Pos/filtered_feature_tables/processed_preferred_Feature_table.tsv  Supernatant_Metabolomics+ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_Hilic_Pos/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Media_18_HILICpos Media_29_HILICpos
#python3 ./Analysis.py ./FeatureTables/Media_Hilic_Neg/filtered_feature_tables/processed_preferred_Feature_table.tsv  Supernatant_Metabolomics- /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_Hilic_Neg/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Media_6_HILICneg

# With group filtering

python3 ./Analysis.py ./FeatureTables/Pellet_RP_Pos/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Pellet_Lipidomics+ /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_RP_Pos/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Pellets_18_Rppos Pellets_9_Rppos
python3 ./Analysis.py ./FeatureTables/Pellet_RP_Neg/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Pellet_Lipidomics- /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_RP_Neg/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Pellets_18_Rpneg Pellets_22_Rpneg
python3 ./Analysis.py ./FeatureTables/Pellet_HILIC_pos/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Pellet_Metabolomics+ /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_HILIC_pos/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Pellets_25_HILICpos Pellets_28_HILICpos
python3 ./Analysis.py ./FeatureTables/Pellet_HILIC_Neg/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Pellet_Metabolomics- /Users/mitchjo/Projects/cyp19a1_targeted_analysis/Sequence_iPSC_Pellets_regularmethods_03252023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Pellet_HILIC_Neg/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Pellets_25_HILICneg

python3 ./Analysis.py ./FeatureTables/Media_RP_Pos/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Supernatant_Lipidomics+ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_RP_Pos/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Media_7_Rppos Media_8_Rppos
python3 ./Analysis.py ./FeatureTables/Media_RP_Neg/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Supernatant_Lipidomics- /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_RP_Neg/annotations/LMSD_annotated_default_empCpds.json ./LMSD.json Media_18_Rpneg Media_7_Rpneg
python3 ./Analysis.py ./FeatureTables/Media_Hilic_Pos/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Supernatant_Metabolomics+ /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_Hilic_Pos/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Media_18_HILICpos Media_29_HILICpos
python3 ./Analysis.py ./FeatureTables/Media_Hilic_Neg/filtered_feature_tables/processed_gfiltering_Feature_table.tsv  Supernatant_Metabolomics- /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/Sequence_iPSC_media_regularmethods_03182023.csv /Users/mitchjo/Projects/MorPHiC_cyp19a1_prlr_ptn/FeatureTables/Media_Hilic_Neg/annotations/HMDB_annotated_default_empCpds.json ./HMDB.json Media_6_HILICneg
