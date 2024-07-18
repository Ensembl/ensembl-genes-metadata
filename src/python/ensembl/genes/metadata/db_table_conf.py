# db_table_conf.py

TABLE_CONF = {
    "assembly": {"method": "per_col", "dkey":"None", "ukey": "None"},
    "organism": {"method": "per_col", "dkey":"None", "ukey": "None"},
    "species": {"method": "per_col", "dkey":"None", "ukey": "None"},
    "assembly_metrics" :  {"method": "per_row", "dkey":"assembly_id", "ukey": "None"},
    "bioproject_lineage": {"method": "per_row_key", "dkey":"assembly_id", "ukey": "None"},
    }