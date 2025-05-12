import subprocess
import os


def extract_largest_file(zip_file, file_pattern, species_name, file_type, download_dir):
    # List all files in the zip archive that match the pattern
    list_command = f"unzip -l {zip_file} {file_pattern}"
    list_output = subprocess.run(list_command, shell=True, text=True, capture_output=True, check=True)
    lines = list_output.stdout.split('\n')

    # Parse the list output to find the largest file
    max_size = 0
    largest_file = None
    for line in lines:
        parts = line.strip().split()
        if len(parts) > 3 and parts[-1].endswith('.fna'):
            size = int(parts[0])
            filename = parts[-1]
            if size > max_size:
                max_size = size
                largest_file = filename

    # Extract the largest file from the zip archive and rename it
    if largest_file:
        new_filename = f"{species_name.replace(' ', '_')}_{file_type}.fna"
        extract_command = f"unzip -p {zip_file} {largest_file} > {os.path.join(download_dir, new_filename)}"
        subprocess.run(extract_command, shell=True, check=True)
        return os.path.join(download_dir, new_filename)
    else:
        print(f"No file matched the pattern {file_pattern} in {zip_file}")
        return None


def download_data(species_name, download_dir, mode):
    # Define filenames for output

    if mode == "orthology":
        rna_zip_file = f"{species_name.replace(' ', '_')}_rna.zip"
        cds_zip_file = f"{species_name.replace(' ', '_')}_cds.zip"
        fasta_filename =  "user_species.fasta"
        user_orthofinder_dir = download_dir + "/user_orthofinder"  # Target directory for the final renamed file

        # Commands to download RNA and CDS
        rna_command = [
            "datasets", "download", "genome", "taxon", species_name, "--include", "rna",
            "--filename", os.path.join(download_dir, rna_zip_file)
        ]
        cds_command = [
            "datasets", "download", "genome", "taxon", species_name, "--include", "cds",
            "--filename", os.path.join(download_dir, cds_zip_file)
        ]

        # Try downloading RNA first, then CDS if RNA fails
        try:
            print(f"Attempting to download RNA for {species_name}...")
            subprocess.run(rna_command, text=True, capture_output=True, check=True)
            print(f"RNA download successful for {species_name}. Extracting...")

            # Attempt extraction and log the extracted file path
            extracted_file = extract_largest_file(os.path.join(download_dir, rna_zip_file), '*rna.fna*', species_name,
                                                  'rna', download_dir)
            print(f"Extracted file path returned: {extracted_file}")

            # Renaming extracted file if it exists
            if extracted_file and os.path.exists(extracted_file):
                os.rename(extracted_file, os.path.join(user_orthofinder_dir, fasta_filename))
                print(f"File renamed and moved to {os.path.join(user_orthofinder_dir, fasta_filename)}")

                # Cleanup: remove RNA zip file
                os.remove(os.path.join(download_dir, rna_zip_file))
                print(f"Removed {rna_zip_file}")

            else:
                print("Extraction reported success, but no file was found for renaming.")

        except subprocess.CalledProcessError as e:
            print(f"RNA download failed for {species_name}. Error: {e}")

            try:
                print(f"Attempting to download CDS for {species_name}...")
                subprocess.run(cds_command, text=True, capture_output=True, check=True)
                print(f"CDS download successful for {species_name}. Extracting...")

                # Attempt extraction and log the extracted file path
                extracted_file = extract_largest_file(os.path.join(download_dir, cds_zip_file), '*cds_from_genomic.fna*',
                                                      species_name, 'cds', download_dir)
                print(f"Extracted file path returned: {extracted_file}")

                # Renaming extracted file if it exists
                if extracted_file and os.path.exists(extracted_file):
                    os.rename(extracted_file, os.path.join(user_orthofinder_dir, fasta_filename))
                    print(f"File renamed and moved to {os.path.join(user_orthofinder_dir, fasta_filename)}")

                    # Cleanup: remove CDS zip file
                    os.remove(os.path.join(download_dir, cds_zip_file))
                    print(f"Removed {cds_zip_file}")

                else:
                    print("Extraction reported success, but no file was found for renaming.")

            except subprocess.CalledProcessError as e:
                print(f"Failed to download CDS for {species_name} as well. Error: {e}")

    #### OFF-target

    elif mode == "off_target":

        rna_zip_file = f"{species_name.replace(' ', '_')}_rna.zip"

        cds_zip_file = f"{species_name.replace(' ', '_')}_cds.zip"

        fasta_filename = "user_index.fasta"

        user_orthofinder_dir = download_dir + "/efficiency"  # Target directory for the final renamed file

        # Commands to download RNA and CDS

        rna_command = [

            "datasets", "download", "genome", "taxon", species_name, "--include", "rna",

            "--filename", os.path.join(download_dir, rna_zip_file)

        ]

        cds_command = [

            "datasets", "download", "genome", "taxon", species_name, "--include", "cds",

            "--filename", os.path.join(download_dir, cds_zip_file)

        ]

        # Try downloading RNA first, then CDS if RNA fails

        try:

            print(f"Attempting to download RNA for {species_name}...")

            subprocess.run(rna_command, text=True, capture_output=True, check=True)

            print(f"RNA download successful for {species_name}. Extracting...")

            # Attempt extraction and log the extracted file path

            extracted_file = extract_largest_file(os.path.join(download_dir, rna_zip_file), '*rna.fna*', species_name,

                                                  'rna', download_dir)

            print(f"Extracted file path returned: {extracted_file}")

            # Renaming extracted file if it exists

            if extracted_file and os.path.exists(extracted_file):

                os.rename(extracted_file, os.path.join(user_orthofinder_dir, fasta_filename))

                print(f"File renamed and moved to {os.path.join(user_orthofinder_dir, fasta_filename)}")

                # Cleanup: remove RNA zip file

                os.remove(os.path.join(download_dir, rna_zip_file))

                print(f"Removed {rna_zip_file}")


            else:

                print("Extraction reported success, but no file was found for renaming.")


        except subprocess.CalledProcessError as e:

            print(f"RNA download failed for {species_name}. Error: {e}")

            try:

                print(f"Attempting to download CDS for {species_name}...")

                subprocess.run(cds_command, text=True, capture_output=True, check=True)

                print(f"CDS download successful for {species_name}. Extracting...")

                # Attempt extraction and log the extracted file path

                extracted_file = extract_largest_file(os.path.join(download_dir, cds_zip_file),
                                                      '*cds_from_genomic.fna*',

                                                      species_name, 'cds', download_dir)

                print(f"Extracted file path returned: {extracted_file}")

                # Renaming extracted file if it exists

                if extracted_file and os.path.exists(extracted_file):

                    os.rename(extracted_file, os.path.join(user_orthofinder_dir, fasta_filename))

                    print(f"File renamed and moved to {os.path.join(user_orthofinder_dir, fasta_filename)}")

                    # Cleanup: remove CDS zip file

                    os.remove(os.path.join(download_dir, cds_zip_file))

                    print(f"Removed {cds_zip_file}")


                else:

                    print("Extraction reported success, but no file was found for renaming.")


            except subprocess.CalledProcessError as e:

                print(f"Failed to download CDS for {species_name} as well. Error: {e}")


def ncbi_new_species_all(user_dir, species_name):
    # Create the required directories
    os.makedirs(user_dir , exist_ok=True)

    os.makedirs((user_dir + "/user_orthofinder" ), exist_ok=True)

    os.makedirs((user_dir + "/user_orthofinder"+ "/orthofinder_input"), exist_ok=True)

    mode = "orthology"


    # Set directory paths


    # Download and process data
    download_data(species_name, user_dir, mode)


def ncbi_new_species_off_target_all(user_dir, species_name):
    # Create the required directories
    os.makedirs(user_dir, exist_ok=True)

    mode = "off_target"

    # Set directory paths
    species_name = species_name.replace("_", " ")
    # Download and process data
    download_data(species_name, user_dir, mode)
# Example usage
#ncbi_new_species_all("test", "Psylliodes chrysocephala")
