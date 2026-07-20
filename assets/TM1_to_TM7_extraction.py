from Bio.Seq import Seq                     # Import Seq class to create biological sequence objects
from Bio.SeqRecord import SeqRecord         # Import SeqRecord class to store sequence + ID + description
from Bio import SeqIO                       # Import SeqIO to read and write FASTA files


# =====================================================
# 1️⃣ FUNCTION TO LOAD FASTA INTO A DICTIONARY
# =====================================================

def load_fasta(fasta_file):                                   # Define a function that reads a FASTA file

    sequence_dict = {}                                         # Create empty dictionary to store sequences

    print("Reading FASTA file...")                             # Print message to show progress

    for record in SeqIO.parse(fasta_file, "fasta"):            # Loop through each sequence in FASTA file

        accession = record.id                                  # Get FASTA header (accession ID)

        print(f"Loaded sequence: {accession}")                 # Print which sequence was loaded

        sequence_dict[accession] = record                       # Store sequence in dictionary using accession as key

    print(f"Total sequences loaded: {len(sequence_dict)}\n")   # Print total number of sequences

    return sequence_dict                                        # Return dictionary to main program


# =====================================================
# 2️⃣ FUNCTION TO PARSE GFF3 FILE
# =====================================================

def parse_gff3(gff_file):                                      # Define function to read GFF3 file

    tm_dict = {}                                               # Create empty dictionary for TM coordinates

    print("Reading GFF3 file...")                              # Print progress message

    with open(gff_file, "r") as file:                          # Open GFF3 file in read mode

        for line in file:                                      # Loop through each line in file

            line = line.strip()                                 # Remove newline and surrounding spaces

            if not line or line.startswith("#") or line.startswith("//"):
                continue                                        # Skip empty lines and comment lines

            parts = line.split("\t")                            # Split line into columns using TAB

            if len(parts) < 4:                                  # Make sure there are at least 4 columns
                continue                                        # Skip invalid lines

            accession = parts[0]                                # Column 1 = accession ID
            feature_type = parts[1]                             # Column 2 = feature type (TMhelix, inside, outside)
            start = int(parts[2])                               # Column 3 = start position (convert to integer)
            end = int(parts[3])                                 # Column 4 = end position (convert to integer)

            if feature_type == "TMhelix":                       # Only process TMhelix rows

                print(f"Found TMhelix: {accession} {start}-{end}")  # Print detected TM region

                if accession not in tm_dict:                    # If accession not yet stored
                    tm_dict[accession] = []                     # Create empty list for this accession

                tm_dict[accession].append((start, end))         # Add (start, end) tuple to list

    print(f"Total sequences with TM annotations: {len(tm_dict)}\n")  # Print total annotated sequences

    return tm_dict                                              # Return TM dictionary


# =====================================================
# 3️⃣ FUNCTION TO EXTRACT AND CONCATENATE TM1–TM7
# =====================================================

def extract_tm_regions(sequence_dict, tm_dict):                # Define function to extract TM regions

    output_records = []                                         # Create empty list to store new TM-only sequences

    print("Extracting TM regions...\n")                        # Print progress message

    for accession, tm_list in tm_dict.items():                  # Loop through each accession in TM dictionary

        print(f"Processing: {accession}")                       # Print which accession is being processed

        if accession not in sequence_dict:                      # Check if accession exists in FASTA
            print("⚠ Not found in FASTA — skipping\n")         # Warning if not found
            continue                                            # Skip to next accession

        full_sequence = str(sequence_dict[accession].seq)       # Get full protein sequence as string

        concatenated_tm = ""                                    # Create empty string to store TM1–TM7

        tm_list_sorted = sorted(tm_list, key=lambda x: x[0])    # Sort TM regions by start coordinate

        print(f"TM regions: {tm_list_sorted}")                  # Print sorted TM list

        for i, (start, end) in enumerate(tm_list_sorted[:7]):   # Loop through first 7 TM helices only

            print(f"Extracting TM{i+1}: {start}-{end}")         # Print which TM helix is extracted

            start_index = start - 1                              # Convert 1-based GFF start to 0-based Python index
            end_index = end                                      # End index remains same in slicing

            if 0 <= start_index < len(full_sequence) and 0 < end_index <= len(full_sequence):
                segment = full_sequence[start_index:end_index]   # Extract TM segment from sequence
                concatenated_tm += segment                       # Add segment to concatenated string
            else:
                print("⚠ Coordinates out of bounds — skipping this TM")  # Warning if coordinates invalid

        if concatenated_tm:                                      # If at least one TM region extracted

            new_record = SeqRecord(                              # Create new SeqRecord object
                Seq(concatenated_tm),                             # Create new sequence object from TM string
                id=f"{accession}_TM1-7",                          # New FASTA header
                description=""        # Description field
            )

            output_records.append(new_record)                     # Add new record to output list

        print()                                                  # Print blank line for readability

    print(f"Total TM-concatenated sequences: {len(output_records)}\n")  # Print final count

    return output_records                                        # Return list of TM-only sequences


# =====================================================
# 4️⃣ MAIN PROGRAM
# =====================================================

fasta_file = "C:/Users/USER/Documents/Master_of_Science_in_Systems_Biology_GRA/Experimental_results/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/renamed_combined_CPB_and_representative_OARs_sequences.fasta"                     # Path to input FASTA file
gff_file = "C:/Users/USER/Documents/Master_of_Science_in_Systems_Biology_GRA/Experimental_results/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/deeptmhmm_results/TMRs.gff3"                  # Path to input GFF3 file
output_file = "C:/Users/USER/Documents/Master_of_Science_in_Systems_Biology_GRA/Experimental_results/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/extracted_TM1-7_output.fasta"                # Path for output FASTA file

sequence_dict = load_fasta(fasta_file)                           # Load FASTA sequences into dictionary

tm_dict = parse_gff3(gff_file)                                   # Parse GFF3 to get TM coordinates

tm_records = extract_tm_regions(sequence_dict, tm_dict)          # Extract TM1–TM7 regions

SeqIO.write(tm_records, output_file, "fasta")                    # Write TM-only sequences to new FASTA file

print("✅ TM extraction completed successfully.")                 # Final completion message