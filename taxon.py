import re
import sys
import argparse
from Bio import Entrez

def get_ncbi_tax(taxon):
  '''Getn NCBI taxonomy'''
  # If the input is a string
  if not re.match(r'\d+', taxon):
    # Get taxonomy ID using Entrez
    taxon2 = '"' + taxon + '"'
    handle = Entrez.esearch(
    db='taxonomy', term=taxon2, rettype='gb', retmode='text')
    record = Entrez.read(handle, validate=False)
    handle.close()
    # If there's no result
    if not record['IdList']:
      sys.exit(
      '[ERROR] The taxon "{}" you provided is invalid. '
      'Please check NCBI Taxonomy'.format(taxon))
      tax_id = record['IdList']
  else:
    tax_id = taxon

  # Now connect NCBI again using the tax_id
  # Entrez.efetch will give you various information
  handle2 = Entrez.efetch(db='taxonomy', id=tax_id, retmode='xml')
  record2 = Entrez.read(handle2, validate=False)
  handle2.close()
  report = ''
  if ('OtherNames' in record2[0]):
    common_name = record2[0]['OtherNames']
    if ('GenbankCommonName' in common_name):
      report = 'common_name:' + common_name['GenbankCommonName'] + ':' + record2[0]['ScientificName'] +'\n'
    elif('CommonName' in common_name):
      report = 'common_name:' + common_name['CommonName'][0] + ':' + record2[0]['ScientificName'] +'\n'
  else:
    report = 'common_name:' + '' + ':' + record2[0]['ScientificName'] +'\n'
  tax_list = record2[0]['LineageEx']
  for tax_element in tax_list:
    report = report + ('{}: {}: {}'.format(
    tax_element['Rank'], tax_element['ScientificName'],  tax_element['TaxId'])) + "\n"
  with open('taxonomy_report.txt', 'a') as writer:
    writer.write(report)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-tid','--taxon_id', help='Taxonomy id to fetch information about', required=True)           
  args = parser.parse_args()
  taxon_id = args.taxon_id
  # Now call the function
  get_ncbi_tax(taxon_id)  # Using tax ID
  #get_ncbi_tax('Lentinula edodes')  # Using tax name
