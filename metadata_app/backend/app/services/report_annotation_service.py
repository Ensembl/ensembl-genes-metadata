import pandas as pd
import logging
from metadata_app.backend.app.services.annotations_service import generate_tables


def generate_report(end_date, start_date, group_name, taxon_id, bioproject_id):
	anno_wide, anno_main = generate_tables(group_name=group_name, taxon_id=taxon_id, bioproject_id=bioproject_id,
	                                       annotation_date=None)
	# Ensure date fields are datetime
	for col in ["last_genebuild_update", "release_date", "date_status_update"]:
		if col in anno_wide.columns:
			anno_wide[col] = pd.to_datetime(anno_wide[col], errors="coerce")

	if end_date:
		end_date = pd.to_datetime(end_date)
		anno_wide = anno_wide[anno_wide["date_status_update"] <= end_date]

	if start_date:
		start_date = pd.to_datetime(start_date)
		anno_wide = anno_wide[anno_wide["date_status_update"] >= start_date]

	# Create tables for charts
	number_of_annotations_raw = (
		anno_wide[['gca', 'gb_status']]
		.groupby('gb_status')
		.size()
		.reset_index(name='count')
	)

	# Transform into desired format
	number_of_annotations = [
		{
			"gb_status": row['gb_status'],
			"count": row['count']
		}
		for _, row in number_of_annotations_raw.iterrows()
	]

	method_report = (
		anno_wide[['gca', 'annotation_method']]
		.groupby('annotation_method')
		.size()
		.reset_index(name='count'))

	num_unique_taxa = anno_wide['lowest_taxon_id'].nunique()
	top_3_taxa = (
		anno_wide.groupby(['scientific_name'])
		.size()
		.reset_index(name='count')
		.sort_values(by='count', ascending=False)
		.head(3)
	)

	project_report = (
		anno_wide[['gca', 'associated_project']]
		.groupby('associated_project')
		.size()
		.reset_index(name='count'))

	if "protein_busco" in anno_wide.columns:
		extracted = anno_wide["protein_busco"].str.extract(r"C:(\d+\.\d+)%")[0]

		# Convert to float, but invalid values → None instead of NaN
		extracted = pd.to_numeric(extracted, errors="coerce")
		# Store cleaned series back in the DataFrame (optional)
		anno_wide["busco_complete"] = extracted
		# Compute average only on valid numbers
		valid = extracted.dropna()

		average_busco = valid.mean() if not valid.empty else "Not available"
	else:
		anno_wide["protein_busco"] = "Not available"
		average_busco = "Not available"

	main_report = anno_wide[
		['associated_project', 'gca', 'genebuilder', 'gb_status', 'ftp', 'latest_annotated', 'protein_busco',
		 'last_genebuild_update', 'release_date']]
	logging.info(f"BUSCo {average_busco}")

	# Transforming out of range float values that are not JSON compliant: nan
	logging.info(f"Transfroming Out of range float values that are not JSON compliant")
	anno_wide = anno_wide.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	method_report = method_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	top_3_taxa = top_3_taxa.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	project_report = project_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)
	main_report = main_report.apply(lambda col: col.fillna("") if col.dtype == "object" else col)

	return anno_wide, number_of_annotations, method_report, num_unique_taxa, top_3_taxa, project_report, average_busco, main_report
