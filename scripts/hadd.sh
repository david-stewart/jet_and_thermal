rm hadd*
hadd hadd_JM.root JM*.root
../../scripts/to_parquet.py hadd_JM.root
