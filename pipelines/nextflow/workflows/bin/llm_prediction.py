#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
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
"""Tissue prediction using LLMs
"""
import time
import argparse
import re
import pandas as pd
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer, pipeline
from huggingface_hub import login
from datasets import Dataset
import pymysql
from tqdm import tqdm


BATCH_SIZE = 50


def connect_to_db(host, user, password, database, port=3306):
    """
    Establishes a connection to the MySQL database.

    :param host: The hostname or IP address of the MySQL server.
    :param user: The username to connect to the database.
    :param password: The password for the user.
    :param db: The name of the database to connect to.
    :param port: The port of the MySQL server (default is 3306).
    :return: A pymysql connection object.
    """
    try:
        connection = pymysql.connect(host=host, user=user, password=password, database=database, port=port)
        return connection
    except pymysql.MySQLError as e:
        print(f"Error connecting to MySQL: {e}")
        return None


def create_prompt(biosample_tissue, sentence) -> str:
    """
    Create a prompt for the LLM model to extract tissue names from the given sentence.
    :param biosample_tissue: The biosample tissue to be used in the prompt.
    :param sentence: The sentence from which to extract tissue names.
    :return: The formatted prompt string.
    """
    return f"""
You are an expert in biomedical text processing specialised in eukaryiotic life ontologies. Extract only valid anatomicial entities (tissues, organs, cells) from biomedical text.
You will be given a biosample tissue or a run description.
If biosample tissue is available, use only that.
If not, then extract tissues from the run description.

Rules:
1. Extract only tissues or body parts (e.g., "liver," "brain," "whole larvae")  and  return them in lower case singular  with no extra text
2. Do NOT extract sequencing or experimental terms (e.g., "Illumina", "NovaSeq", "RNA-seq")
3. STRICTLY IGNORE all taxonomy scientific names and common names of vertebrates, invertebrates, and plants. If an example  matches an animal species only, return "NONE". Example ignored species include (case insensitive): Macaca fascicularis, Macaca nemestrina, Rattus norvegicus, Bos grunniens, Cricetulus griseus, Homo sapiens, Pan troglodytes, Oryctolagus cuniculus, Mus musculus, Gallus gallus, Danio rerio, Drosophila melanogaster, Arabidopsis thaliana, Zea mays, rabbit, sheep, moths, salmon
4. Normalize extracted tissue names: replace underscores (_) with spaces, remove new line, numbers, dates, stopwords, singularize plurals, remove duplicates.
5. If multiple tissues are mentioned, return them in a comma-separated list without repetition.
6. If the sentence contains explicit reference to cancer or disease (e.g., 'glioblastoma', 'tumor sample', 'g2_t2', 'tumor'), add 'CANCER/DISEASE' next to it.
7. Ignore IDs, alphanumeric codes, and random labels (e.g., xz_l_13b, sample_XYZ123)
8  Only output the tissue(s) with no explanation. Do not explain your reasoning. Output only the final result.
9. If no valid tissue is found, return "NONE" with no extra text.

Examples:
- "Illumina NovaSeq 6000 paired end sequencing" -> "NONE"
- "pool of 60 whole larvae" → "whole larvae"
- "plasma" -> "plasma"
- "Illumina NovaSeq paired-end sequencing" → "NONE"
- "liver sample" → "liver"
- "brain and spinal cord" → "brain, spinal cord"
- "RNA extraction from kidney" → "kidney"
- "Illumina HiSeq 4000, Illumina HiSeq 2000 sequencing Unknown_BX820041T0016, GSM8569180r1 Illumina HiSeq 2000 sequencing GSM8569180 gADSCs Control1 Capra hircus RNASeq" → "NONE"
- "nextseq sequencing gsm fomite mesocricetus auratus rnaseq, nextseq sequencing gsm df" -> "NONE"
- "whole_body" -> "whole body"
- "glioblastoma" -> "glioblastoma, CANCER/DISEASE tissue"
- "sars-cov-2, sarscov2 infection" -> "sarscov2, CANCER/DISEASE tissue"
- "devil facial tumour 1" -> "devil facial tumour, CANCER/DISEASE tissue"
- "ecffhchc illumina hiseq sequencing rnaseq acanthochromis polyacanthus juvenile liver illumina hiseq sequencing rnaseq acanthochromis polyacanthus juvenile liver" -> "liver"
- "Illumina HiSeq 2000 paired end sequencing Illumina HiSeq 2000 paired end sequencing of sheep SAMD00115868 gonad DRX121131 Illumina HiSeq 2000 paired end sequencing Illumina HiSeq 2000 paired end sequencing of SAMD00115868" -> "gonad"
- "2cell embryo capra hircus" -> "embryo"
- "34 l 36, 34_l_amy, qbt_l_tea, xz_l_13b,zs_l_iapm_rep,2ramy,mfacad710, macaca fascicularis,equus asinus, ZEA10_2,cyno_5,"tasmanian_devil, sarcophilus harrisii, rnaseq, nan, gsm62, xq mst, xq pu rep,zs bd" -> "NONE"
- "cell line" -> "cell line"
- "pbmcs" -> "pbmc"
- "adult testis, callithrix jacchus" -> "adult testis"  REMOVE SPECIES NAME
- "baboon, spcaf, papio anubis, whart, parietallobeleft" -> "parietal lobe left"
- "rna extraction from kidney" -> "kidney"
- "macaca fascicularis,  Homo sapiens, Pan troglodytes" -> "NONE"
- "the study used rna sequencing to analyze the transcriptome of the whole larvae of the fruit fly drosophila mel" -> "whole larvae"
- "hiseq ten sequencing gsm mfa shc macaca fascicularis rnaseq gsm hiseq ten sequencing" -> "NONE"
- "illumina hiseq sequencing gsm fcg felis catus rnaseq, illumina hiseq sequencing g" -> "NONE"
- "RNASeq" -> "NONE"
- "Callithrix jacchus, NMiPS, RNASeq, Horse F IL1b licensed Equus caballus RNASeq,control DMSO Cricetulus griseus" -> "NONE"
- "adipose 1, adipose1" -> "adipose"
- "heart1 of yak bos grunni" -> "heart"
- "mink testis april 1" -> "mink testis"
- "rattus norvegicus, pan paniscus, macaca nemestrina,Cricetulus griseus, Tibetan sheep" -> "NONE"
- "the longest back muscle 1, muscle17" -> "muscle"
- "Muscle 918 Ovis aries RNASeq" -> "muscle"
- "Illumina HiSeq 2000 sequencing GSM4026140 Mammary gland cells Lv5NC1 Capra hircus RNASeq" -> "mammary gland cells"
- "Illumina HiSeq 2500 sequencing Capra hircus transcriptome 12monthB3" -> "NONE"
- "Illumina NovaSeq 6000 sequencing GSM4297844 singleCell_spheroid_plate2_H6_S251 Pan troglodytes Homo sapiens RNASeq" -> "NONE"
- "rabbit, pig, sheep, chimpanzee, human, mice, rats, baboon, rhesus" ->"NONE"
- "large scale transcriptional analysis of mhc class haplotype diversity in sheep" -> "mhc"
- "root" -> "root"
- "Illumina HiSeq 4000 sequencing GSM6588584 ART_3 Muscle Bos indicus x Bos taurus RNASeq skin" -> "muscle, skin"
- "HiSeq X Ten sequencing GSM4615566 I5332F_Fibroblast Callithrix jacchus RNASeq" -> "Fibroblast"
- "oviduct" -> "oviduct"
- "Illumina NovaSeq 6000 sequencing kidney" -> "kidney"
- "qbt_l_11m, xq_r_ia, 92_l_ial, 34_r_na, qbt_l_g" -> "NONE"
- "CD4 subset of Goat at day 5" -> "CD4"
- "muscle, barnfed goat muscle 3" -> "muscle"
- "midbrainleft" -> "mid brain left"
- "MiscOrgans_START1, Cricetulus griseus, ncRNASeq, whole organism" -> "whole organism"
- "lung10" -> "lung"
- "lowfat_S180_ME_L001_R" -> "low fat"
- "Illumina MiSeq paired end sequencing, NextSeq 500 sequencing" -> "NONE"
- "muscle, CANCER/DISEASE tissue, CANCER/DISEASE tissue, CANCER/DISEASE tissue, CANCER/DISEASE tissue" -> "muscle"
- "subcutaneous adipose 1" -> "subcutaneous adipose"
Based on the above rules and examples, now process the following biosample tissue "{biosample_tissue}". If not exists process run description "{sentence}" 
Answer:
"""


def load_model(token_input: str) -> pipeline:
    """Model loading function.
    This function loads the Mistral-7B-Instruct-v0.3 model from Hugging Face and prepares it for inference.
    It uses the AutoModelForCausalLM and AutoTokenizer classes from the transformers library.
    The model is loaded with gradient checkpointing enabled to save memory.
    The model is set to use half-precision (FP16) for lower memory usage.
    The device is set to GPU if available, otherwise it defaults to CPU.
    The pipeline is created for text generation using the loaded model and tokenizer.
    The function returns the pipeline object for further use.
    This function requires the following libraries:

    Args:
        token_input (str): needed to load the model from Hugging Face.

    Returns:
        pipeline: pipeline object for text generation.
    """
    torch.cuda.empty_cache()  # Clears unused memory before loading
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    model_name = "mistralai/Mistral-7B-Instruct-v0.3"
    tokenizer = AutoTokenizer.from_pretrained(
        model_name, token=token_input
<<<<<<< HEAD
    ) 
=======
    )  # "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"
>>>>>>> a70c544... llm prediction script

    if tokenizer.pad_token_id is None:
        tokenizer.pad_token = tokenizer.eos_token  # Set pad token to EOS
    model = AutoModelForCausalLM.from_pretrained(
        model_name,
        torch_dtype=torch.float16,  # Use FP16 for lower memory
        device_map="auto",  # Distributes layers across available GPUs
    )
    # Enable gradient checkpointing to save memory
    model.gradient_checkpointing_enable()
    pipe = pipeline("text-generation", model=model, tokenizer=tokenizer)  # , device=device)

    return pipe


def define_model_dataset(df, column_input, column_output, token_input: str, batch_size: int) -> pd.DataFrame:
    """Define the model and dataset for LLM prediction.
    This function loads the model and prepares the dataset for inference.
    It uses the Hugging Face pipeline for text generation.
    The function takes a DataFrame, the input column name, and the output column name as arguments.
    It creates a Hugging Face Dataset from the DataFrame and uses the pipeline to generate predictions.
    The function returns the DataFrame with the predictions added to the output column.
    This function requires the following libraries:
    Args:
        df (pd.DataFrame): DataFrame containing the input data.
        column_input (str): Name of the input column in the DataFrame.
        column_output (str): Name of the output column in the DataFrame.
        token_input (str): Needed to load the model from Hugging Face.
    Returns:
        pd.DataFrame: DataFrame with the predictions added to the output column.
    """
    # Load the model
    pipe = load_model(token_input=token_input)
    # Convert Pandas DataFrame to Hugging Face Dataset
    dataset = Dataset.from_pandas(df)
    start_time = time.time()
    responses = pipe(
        dataset[f"{column_input}"],  # Directly pass all prompts
        batch_size=batch_size,
        max_new_tokens=60,
        return_full_text=True,
    )

    # ✅ Assign responses back to DataFrame (keeping order)
    df[f"{column_output}"] = [
        resp[0]["generated_text"].replace(df.iloc[i][f"{column_input}"], "").strip()
        for i, resp in enumerate(responses)
    ]
    print(f"Total time: {time.time() - start_time:.2f} sec")

    return df


def clean_tissue_prediction(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the tissue prediction column in the DataFrame.
    This function processes the tissue prediction column to remove unwanted characters and formats the text.
    It uses regular expressions to clean the text and applies additional formatting rules.
    The function returns the cleaned DataFrame.
    Args:
        df (pd.DataFrame): DataFrame containing the tissue prediction column.
    Returns:
        pd.DataFrame: DataFrame with the cleaned tissue prediction column.
    """
    # df["tissue_prediction"] = df["tissue_prediction"].apply(
    #            lambda x: re.search(r"^(.*)", str(x)).group(1) if re.search(r"^(.*)", str(x)) else x)

    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: re.sub(r"\d+|day", "", str(x)).strip())

    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: re.sub(r"\b[a-zA-Z]\b", "", str(x)).strip("_")
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1) if (m := re.search(r"^(.*?)\s*Based on", str(x))) else x)
    )
    #df["tissue_prediction"] = df["tissue_prediction"].apply(
    #    lambda x: (
    #        re.search(r"^(.*?)\s*Based on", str(x)).group(1) if re.search(r"^(.*?)\s*Based on", str(x)) else x
    #    )
    #)
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1) if (m := re.search(r"^(.*?)\s*Based on the above rules and examples", str(x))) else x)
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1) if (m := re.search(r"^(.*?)\s*If not available:", str(x))) else x)
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1) if (m := re.search(r"^(.*?)\s*If .*? is not available, then:", str(x))) else x)
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1) if (m := re.search(r"^(.*?)\s*->", str(x))) else x)
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1) if (m := re.search(r"^(.*?)\s*If run description is", str(x))) else x)
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: None if "is not a valid tissue" in str(x) else x
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1) if (m := re.search(r"^(.*?)\s*Now process .*? run description", str(x))) else x)
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1).strip() if (m := re.search(r"\(if run description is used\)(.*?)\(if biosample tissue is used\)", str(x))) else x)
    )
    df["tissue_prediction"] = (
        df["tissue_prediction"]
        .astype(str)
        .str.replace('"', "", regex=True)
        .replace("Answer: ", "", regex=True)
        .replace("\n", "", regex=True)
        .replace("Output: ", "", regex=True)
    )
    df.loc[
        df["tissue_prediction"].astype(str).str.contains("NONE", case=False, na=False), "tissue_prediction"
    ] = None
    df.loc[
        df["tissue_prediction"].astype(str).str.contains("nan", case=False, na=False), "tissue_prediction"
    ] = None
    return df


def main() -> None:
    """Module's entry-point."""
    parser = argparse.ArgumentParser(prog="llm_prediction.py", description="Predict tissue using LLMs")
    parser.add_argument(
        "--hugging_face_token",
        type=str,
        required=True,
        help="Hugging Face login token",
    )
    parser.add_argument("--taxonomy_id", default="10116", type=str, required=True, help="Taxonomy ID")
    parser.add_argument("--host",  type=str, required=True, help="Host")
    parser.add_argument("--user",  type=str, required=True, help="User")
    parser.add_argument("--password",  type=str, required=True, help="Password")
    parser.add_argument(
        "--database",
        type=str,
        required=True,
        help="Database",
    )
    parser.add_argument("--port", type=int, required=True, help="Port")
    parser.add_argument("--batch_size", default=32, type=int, required=False, help="Batch size")

    args = parser.parse_args()
    db_config = {
        "host": args.host,
        "user": args.user,
        "password": args.password,
        "database": args.database,
        "port": args.port,
    }
<<<<<<< HEAD
    login(args.hugging_face_token)
=======
    login(args.hugging_face_token)  # "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"
>>>>>>> a70c544... llm prediction script

    # Load the rows to fix (e.g., 150k rows)
    query = (
        "SELECT run_accession,run_description,tissue,sample_tissue  FROM run where taxon_id="
        + args.taxonomy_id
        + " and tissue_prediction is NULL"
    )
    # Connect to the database
    db_connection = connect_to_db(**db_config)
    # Clean the input text
    cursor = db_connection.cursor()
    if db_connection:
        df = pd.read_sql(query, db_connection)
    print(f"Loaded {len(df)} rows.")

    start_time = time.time()
    df["llm_prompt"] = ""
    df["tissue_prediction"] = ""
    print(df.head())
    for i in range(0, len(df)):
        input_sentence = " ".join(df.iloc[i][["tissue", "run_description"]].dropna().values.astype(str))
        prompt = create_prompt(df["sample_tissue"][i], input_sentence)
        df.at[i, "llm_prompt"] = prompt
    define_model_dataset(df, "llm_prompt", "tissue_prediction", args.hugging_face_token, args.batch_size)
    print(f"Total time: {time.time() - start_time:.2f} sec")
    df = clean_tissue_prediction(df)
    #df.to_csv(
    #    "/nfs/production/flicek/ensembl/genebuild/ftricomi/script/learning-francesca/data/minitest_cleaned.csv",
    #    sep="\t",
    #) 

    # Update back to DB
    update_query = "UPDATE run set tissue_prediction = %s where run_accession = %s"
    counter = 0
    for _, row in tqdm(df.iterrows(), total=len(df)):
        try:
            run_acc = row["run_accession"]
            tissue = row["tissue_prediction"]
            print(f"Updating row {row['run_accession']}")
            cursor.execute(update_query, (tissue, run_acc))
            counter += 1
            # Commit every BATCH_SIZE rows
            if counter % BATCH_SIZE == 0:
                db_connection.commit()
                print(f"Committed {counter} rows...")
        except Exception as e:
            print(f"Failed to update row {row['run_accession']}: {e}")

    # Commit changes
    db_connection.commit()
    cursor.close()
    db_connection.close()
    print("✅ Update complete.")


if __name__ == "__main__":
    main()
