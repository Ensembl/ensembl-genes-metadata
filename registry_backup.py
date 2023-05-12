""" Script to backup the meta database """

import re
import os
import datetime
import pymysql as MySQLdb
import argparse

def copy_db(backup_file,host,db,user,port):
    con = MySQLdb.connect(host=host, db=db, user=user, port=int(port), passwd='')
    cur = con.cursor()

    cur.execute("SHOW TABLES")
    data = ""
    tables = []
    for table in cur.fetchall():
        tables.append(table[0])
    for table in tables:
        data += "DROP TABLE IF EXISTS `" + str(table) + "`;"

        cur.execute("SHOW CREATE TABLE `" + str(table) + "`;")
        data += "\n" + str(cur.fetchone()[1]) + ";\n\n"

        cur.execute("SELECT * FROM `" + str(table) + "`;")
        for row in cur.fetchall():
            data += "INSERT INTO `" + str(table) + "` VALUES("
            first = True
            for field in row:
                if not first:
                    data += ', '
                if field is None:
                    value = 'NULL'
                    data +=  str(value)
                elif field == 0:
                    data += "'" + str(field) + "'"
                elif not field:
                    value = 'NULL'
                    data +=  str(value)
                else:
                    field = re.sub("[$@&'\"#~^!?]","",str(field)) 
                    data += "'" + str(field) + "'"
                first = False


            data += ");\n"
        data += "\n\n"

    FILE = open(backup_file,"w")
    FILE.writelines(data)
    FILE.close()


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-p','--port', help='Port number for host', required=True)
  parser.add_argument('-u','--user', help='Mysql user', required=True)
  parser.add_argument('-db','--dbname', help='Database to be backed up', required=True)
  parser.add_argument('-host','--server', help='Host server for database', required=True)
  parser.add_argument('-bkup','--backup_file', help='File name to hold db backup', required=True)
  args = parser.parse_args()
  bkup_file = args.backup_file
  server = args.server
  db = args.dbname
  user = args.user
  port = args.port
  port.strip
  copy_db(bkup_file,server,db,user,port)
