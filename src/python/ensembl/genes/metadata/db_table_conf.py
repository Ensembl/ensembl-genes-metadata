# db_table_conf.py

TABLE_CONF = {
    "assembly": {"method": "per_col", "dkey":"None", "ukey": "None"},
    "organism": {"method": "per_col", "dkey":"assembly_id", "ukey": "organism_id"},
    "species": {"method": "per_col", "dkey":"None", "ukey": "None"},
    "assembly_metrics" :  {"method": "per_row", "dkey":"assembly_id", "ukey": "None"},
    "bioproject": {"method": "per_row_key", "dkey":"assembly_id", "ukey": "None"},
    "group_assembly": {"method": "per_col", "dkey":"assembly_id", "ukey": "assembly_id"},
    }