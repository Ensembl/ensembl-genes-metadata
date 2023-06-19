# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

import argparse
import os
import requests

def slack_reporting(report):
    payload="{\"channel\": \"@denye\", \"username\": \"registry_messenger\", \"text\": \"" + report  +"\"}"
    url = os.getenv('slack_token')
    headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
    r = requests.post(url, data=payload, headers=headers)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--msg', help='Message to report', required=True)
  args = parser.parse_args()
  msg = args.msg
  slack_reporting(msg)
