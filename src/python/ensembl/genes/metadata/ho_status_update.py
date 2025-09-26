import argparse
import pandas as pd
import pymysql
import logging

# Configure logging
logging.basicConfig(
	level=logging.INFO,
	format="%(asctime)s [%(levelname)s] %(message)s",
	handlers=[
		logging.FileHandler("genebuild_status_updater.log"),
		logging.StreamHandler()
	]
)
logger = logging.getLogger(__name__)


def mysql_fetch_data(query, host, user, password, port, database):
	"""
	Execute a MySQL query and return the results as a list of dictionaries.

	Args:
		query (str): SQL query to execute.
		host (str): MySQL host.
		user (str): MySQL username.
		password (str): MySQL password.
		port (int): MySQL port.
		database (str): Database name.

	Returns:
		list[dict]: List of rows as dictionaries.
	"""
	try:
		conn = pymysql.connect(
			host=host,
			user=user,
			password=password,
			port=port,
			database=database,
			cursorclass=pymysql.cursors.DictCursor
		)

		with conn:
			with conn.cursor() as cursor:
				cursor.execute(query)
				results = cursor.fetchall()
				logger.info(f"Query returned {len(results)} rows.")
				return results

	except pymysql.MySQLError as e:
		logger.error(f"MySQL error: {e}")
		return []


def get_genebuild_status():
	"""
	Fetch genebuild status from the registry excluding 'archive' and 'live'.


	Returns:
		dict or None: First row of results as a dictionary, or None if no results.
	"""
	try:
		registry_query = """
            SELECT 
                gca_accession AS accession, 
                genebuild_status_id,
                gb_status
            FROM genebuild_status 
			WHERE gb_status NOT IN ('archive', 'live')
        """

		gb_status = mysql_fetch_data(
			registry_query,
			host="mysql-ens-genebuild-prod-1",
			user="ensro",
			port=4527,
			database="gb_assembly_metadata",
			password=""
		)

		logger.info(f"Found {len(gb_status)} entries in genebuild_status table.")
		gb_status = pd.DataFrame(gb_status)
		return gb_status

	except pymysql.Error as err:
		logger.error("MySQL error: %s", err)
		return {}


def check_status_production_db(gca_tuple):
	try:
		production_query = f"""
            SELECT 
                assembly.accession,
                dataset.status,
                ensembl_release.release_date
            FROM assembly
            JOIN genome
                ON assembly.assembly_id = genome.assembly_id
            JOIN genome_dataset
                ON genome.genome_id = genome_dataset.genome_id
            JOIN dataset
                ON genome_dataset.dataset_id = dataset.dataset_id
            JOIN dataset_attribute
                ON dataset.dataset_id = dataset_attribute.dataset_id
            JOIN dataset_source
                ON dataset.dataset_source_id = dataset_source.dataset_source_id
            JOIN ensembl_release on genome_dataset.release_id = ensembl_release.release_id
            WHERE dataset.name = "genebuild"
                AND assembly.accession IN {gca_tuple}
              AND genome_dataset.is_current = 1
        """

		production_status = mysql_fetch_data(
			production_query,
			host="mysql-ens-production-1",
			user="ensro",
			port=4721,
			database="ensembl_genome_metadata",
			password=""
		)
		production_status = pd.DataFrame(production_status)
		production_status = production_status.drop_duplicates(subset='accession', keep='first')
		print(production_status)

		logger.info(f"Found {len(production_status)} entries in production table.")
		return production_status

	except pymysql.Error as err:
		logger.error("MySQL error: %s", err)
		return {}


def get_status_updates(merged_df):
	"""
	Determine which GCAs need their genebuild status updated.

	Args:
		merged_df (pd.DataFrame): Merged DataFrame with columns:
			- 'genebuild_status_id'
			- 'gb_status' (current registry status)
			- 'status' (production status)
			- 'release_date'

	Returns:
		pd.DataFrame: DataFrame with columns:
			- 'genebuild_status_id'
			- 'gb_status_new'
			- 'release_date'
		Only rows where status changes are returned.
	"""

	# Copy to avoid modifying the original
	df = merged_df.copy()

	# Initialize a new column for updated status
	df['gb_status_new'] = df['gb_status']

	# Condition 1: Skip Faulty
	condition_faulty = df['status'] == 'Faulty'

	# Condition 2: Released in production but not 'live' in registry
	condition_released = (df['status'] == 'Released') & (df['gb_status'] != 'live')
	df.loc[condition_released, 'gb_status_new'] = 'live'

	# Copy release_date for released entries
	df.loc[condition_released, 'release_date_new'] = df.loc[condition_released, 'release_date']

	# Condition 3: Processing, Processed, Submitted but not released
	condition_handed_over = df['status'].isin(['Processed', 'Processing', 'Submitted']) & (
				df['gb_status'] != 'handed_over')
	df.loc[condition_handed_over, 'gb_status_new'] = 'handed_over'

	# Only keep rows where status actually changes
	updated_df = df[df['gb_status'] != df['gb_status_new']].copy()

	# Keep only relevant columns
	updated_df = updated_df[['genebuild_status_id', 'gb_status_new', 'release_date_new']]
	logging.info(f"Found {len(updated_df)} annotations to be updated.")

	return updated_df


def update_genebuild_status(updated_df, password):
	"""
	Update the genebuild_status table with new statuses and release dates.

	Args:
		updated_df (pd.DataFrame): DataFrame with columns:
			- genebuild_status_id
			- gb_status_new
			- release_date_new (optional, only for 'live')
		password: MySQL connection info
	"""
	connection = pymysql.connect(
		host="mysql-ens-genebuild-prod-1",
		user="enadmin",
		password=password,
		port=4527,
		database="gb_assembly_metadata",
		autocommit=True
	)

	try:
		with connection.cursor() as cursor:
			for _, row in updated_df.iterrows():
				status = row['gb_status_new']
				genebuild_status_id = row['genebuild_status_id']
				release_date = row.get('release_date_new', None)

				if status == 'live' and pd.notnull(release_date):
					sql = """
                        UPDATE genebuild_status
                        SET gb_status = %s,
                            release_date = %s
                        WHERE genebuild_status_id = %s
                    """
					cursor.execute(sql, (status, release_date, genebuild_status_id))
				else:
					sql = """
                        UPDATE genebuild_status
                        SET gb_status = %s
                        WHERE genebuild_status_id = %s
                    """
					cursor.execute(sql, (status, genebuild_status_id))

		logging.info(f"Updated {len(updated_df)} genebuild_status rows.")

	except pymysql.Error as e:
		logging.error("MySQL error: %s", e)
		raise

	finally:
		connection.close()


def main(password):
	logger.info("Fetching genebuild status from registry.")
	gb_status = get_genebuild_status()

	# Make sure we have an 'accession' column
	if 'accession' not in gb_status.columns:
		logger.error("No 'accession' column found in registry data.")
		return

	# Convert GCAs to tuple for SQL IN clause
	gca_tuple = tuple(gb_status['accession'].unique())
	if len(gca_tuple) == 1:
		gca_tuple = (gca_tuple[0],)

	logger.info(f"Getting GCA status from production DB.")
	production_status = check_status_production_db(gca_tuple)

	# Ensure consistent column names
	if 'accession' not in production_status.columns:
		logger.error("Expected columns missing from production DB query.")
		return

	# Merge the two dataframes on 'accession'
	merged_df = pd.merge(
		gb_status,
		production_status,
		on='accession',
		how='inner',
		suffixes=('_registry', '_production')
	)

	updates = get_status_updates(merged_df)

	update_genebuild_status(updates, password)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Update genebuild_status table from production DB.")
	parser.add_argument("-p", "--password", required=True, help="MySQL password for write user")
	args = parser.parse_args()

	main(args.password)





