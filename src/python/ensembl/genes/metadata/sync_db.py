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
from pathlib import Path

""" Method to restore a copy of the updated database across alternate servers """
def restore_db(backup_file,host,user,port,dbname,password):
    cfg_path = Path(backup_file)
    # Backup updated database for restore to fail safe servers
    host = host.split('.')
    mysql_dump = host[0] + ' mysqldump --databases '+ db + ' > ' + backup_file
    os.system(mysql_dump)

    #make sure backup file exists before attempting to process it
    print('Backup file is '+str(cfg_path))
    if cfg_path.is_file():
        sync_server1 = 'mysql -h '+ os.getenv('GBS5') +' -P ' + os.getenv('GBP5') + ' -u' + user +' -p' +password +'  '+ db + ' < ' + str(cfg_path)
        os.system(sync_server1)
        sync_server2 = 'mysql -h '+ os.getenv('GBS6') +' -P ' + os.getenv('GBP6') + ' -u' + user +' -p' +password +'  '+ db + ' < ' + str(cfg_path)
        os.system(sync_server2)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--p', help='Password for host', required=True)
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
  password = args.p
  restore_db(bkup_file,server,user,port,db,password)
