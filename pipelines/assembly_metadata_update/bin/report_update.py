#!/usr/bin/env python3

#  See the NOTICE file distributed with this work for additional information
#  regarding copyright ownership.
#
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import logging
import argparse
import json
import pymysql
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError


def execute_query(query, db_params):
    conn = pymysql.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    return result

def fetching_user(accession, metadata_params, slack_users):

    query_gb_user = f"""SELECT genebuilder, gb_status FROM genebuild_status
                            WHERE gca_accession = '{accession}' and 
                            gb_status in ('in_progress', 'check_busco', 'completed', 'pre_released')"""
    try:
        genebuilder, gb_status = execute_query(query_gb_user, metadata_params)[0]
        slack_id = slack_users[genebuilder]

    except:
        logging.info(f"No available user for {accession}")
        slack_id = None
        gb_status = None

    return slack_id, gb_status

def slack_message(accession, gb_status, check_type, previous_value, current_value):

    custom_message = f"""Accession {accession} with status: {gb_status} has been updated.
    {check_type}: {previous_value} -> {current_value}
    """
    logging.info(f"Slack message: {custom_message}")

    return custom_message

def slack_communication_dm(slack_bot_token, slack_id, custom_message):

    client = WebClient(token=slack_bot_token)

    message = f"""
    Hey <@{slack_id}>. 
    {custom_message}
    """

    resp = client.conversations_open(users=slack_id)
    dm_channel_id = resp["channel"]["id"]
    try:
        client.chat_postMessage(
            channel=dm_channel_id,
            text=message,
            link_names=True, 
        )
    
    except SlackApiError as e:
        logging.error(f"Slack error: {e.response['error']}")


def main():
    """ Module's entry point
    """

    logging.basicConfig(filename="slack_reporting.log", level=logging.DEBUG, filemode='w',
                    format="%(asctime)s:%(levelname)s:%(message)s")
    
    parser = argparse.ArgumentParser(prog='report_update.py',
                                    description="Retrieve metadata from NCBI API for a given GCA accession and store it in JSON files to be inserted in the database.")
    
    parser.add_argument('--accession',
                        type=str,
                        required=True,
                        help='GCA accession to retrieve metadata.')
    parser.add_argument('--check_type',
                        type=str,
                        required=True,
                        help='Type of update check. Possible values: asm_status, refseq, etc.')
    parser.add_argument('--previous_value',
                        type=str,
                        required=True,
                        help='Original value previously stored in the registry.')
    parser.add_argument('--current_value',
                        type=str,
                        required=True,
                        help='Current value previously stored in the registry.')
    parser.add_argument('--slack_params',
                        type=json.loads,
                        required=True,
                        help='Slack bot connection params.')
    parser.add_argument('--slack_users',
                        type=str,
                        required=True,
                        help='Path to json file with username and slack IDs.')
    parser.add_argument('--metadata_params',
                        type=json.loads,
                        required=True,
                        help='JSON/Dict format of gb_assembly_metadata connections params.')
    

    args = parser.parse_args()
    logging.info(args)

    with open(args.slack_users) as json_file:
        slack_users = json.load(json_file)
    
    slack_id, gb_status = fetching_user(args.accession, args.metadata_params, slack_users)

    if slack_id is not None and gb_status is not None:
    
        custom_message = slack_message(args.accession, gb_status, args.check_type, args.previous_value, args.current_value)

        slack_communication_dm(args.slack_params['slack_bot_token'], slack_id, custom_message)
        
        logging.info(custom_message)


if __name__ == '__main__':
    main()