#!/usr/bin/env nextflow
/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

// module description 
process ADAPTER_TRIMMING {
    scratch true
    label 'default'
    tag 'trimming '

    input:
    tuple val(taxon_id), val(run_accession), val(pairedFastqFiles)
    set pair1, pair2 from pairedFastqFiles

    output:
    tuple (taxon_id, run_accession, trimmedFastqFiles) into runTrimmedFastqs
    //no need to emit the path because the subsampled files will be in the run_accession dir _1_1

    script:
    """
    // Construct the command based on whether last_date is provided
    def pythonScript = file("$projectDir/src/python/ensembl/genes/ ")
    def command = "python ${pythonScript} "

    // Execute the Python script
    def process = command.execute()
    process.waitFor()

    // Check if the script execution was successful
    if (process.exitValue() != 0) {
        throw new RuntimeException("Error executing Python script: ${pythonScript}")
    }

    // Read the output of the Python script
    def output = []
    process.inputStream.eachLine { line ->
            output.add(line.trim())
        }

    // Define the sampledFastqFiles array
    def trimmedFastqFiles = []

    // Add the captured paths from the output array to sampledFastqFiles
    output.each { path ->
        trimmedFastqFiles.add(new File(path))
    }

    // Check if the correct number of files were captured
    if (trimmedFastqFiles.size() != 2) {
        throw new RuntimeException("Expected two paths from Python script, but received ${sampledFastqFiles.size()}")
    }

    // Emit the tuple
    emit(taxon_id, run_accession, trimmedFastqFiles)
    """
}

/*
def run_trimming(
    output_dir: Path,
    short_read_fastq_dir: Path,
    delete_pre_trim_fastq: bool = False,
    num_threads: int = 1,
    trim_galore_bin="trim_galore",
) -> None:
    """
    Trim list of short read fastq files.
    Args:
        output_dir : Working directory path.
        short_read_fastq_dir : Short read directory path.
        delete_pre_trim_fastq : Removing original fastq file post trimming. Defaults to False.
        num_threads : Number of threads.
        trim_galore_bin : Software path.
    """
    check_exe(trim_galore_bin)
    trim_dir = create_dir(output_dir, "trim_galore_output")

    fastq_file_list = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    fastq_file_list = [
        path for file_type in file_types for path in Path(short_read_fastq_dir).rglob(file_type)
    ]
    fastq_file_list = _create_paired_paths(fastq_file_list)

    trim_galore_cmd = [
        str(trim_galore_bin),
        "--illumina",
        "--quality",
        "20",
        "--length",
        "50",
        "--output_dir",
        str(trim_dir),
    ]

    pool = multiprocessing.Pool(int(num_threads))  # pylint:disable=consider-using-with
    for fastq_file in fastq_file_list:
        if delete_pre_trim_fastq:
            fastq_file.unlink()
        pool.apply_async(
            multiprocess_trim_galore,
            args=(
                trim_galore_cmd,
                fastq_file,
                trim_dir,
            ),
        )

    pool.close()
    pool.join()

    trimmed_fastq_list = trim_dir.glob("*.fq.gz")

    for trimmed_fastq_path in trimmed_fastq_list:
        logging.info("Trimmed file path: %s", str(trimmed_fastq_path))
        sub_patterns = re.compile(r"|".join(("_val_1.fq", "_val_2.fq", "_trimmed.fq")))
        updated_file_path_name = sub_patterns.sub(".fq", trimmed_fastq_path.name)
        updated_file_path = short_read_fastq_dir / updated_file_path_name
        logging.info("Updated file path: %s", str(updated_file_path))
        trimmed_fastq_path.rename(updated_file_path)

    files_to_delete_list: List[Path] = []
    for file_type in file_types:
        files_to_delete_list.extend(short_read_fastq_dir.glob(file_type))

    for file_to_delete in files_to_delete_list:
        file_to_delete.unlink()
        */