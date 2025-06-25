# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" Script to backup the meta database """

import os
import argparse

def copy_db(backup_file,host,db,user,port):
    host = host.split('.')
    mysql_dump = host[0] + ' mysqldump --databases '+ db + ' | gzip > ' + backup_file + '.gz' 
    os.system(mysql_dump)
    

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--port', help='Port number for host', required=True)
  parser.add_argument('--user', help='Mysql user', required=True)
  parser.add_argument('--dbname', help='Database to be backed up', required=True)
  parser.add_argument('--server', help='Host server for database', required=True)
  parser.add_argument('--backup_file', help='File name to hold db backup', required=True)
  args = parser.parse_args()
  bkup_file = args.backup_file
  server = args.server
  db = args.dbname
  user = args.user
  port = args.port
  port.strip
  copy_db(bkup_file,server,db,user,port)
