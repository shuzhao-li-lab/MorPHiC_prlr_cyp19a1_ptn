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
