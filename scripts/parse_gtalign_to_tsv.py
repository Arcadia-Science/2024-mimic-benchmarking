import argparse
import csv
import re

def parse_gtalign_output(input_file, output_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Parse the header and data sections
    header_line = [
        "num", "reference_description", "tmscore_reference_length_normalized", "tmscore_query_length_normalized", "rmsd", "aligned_residues", "query_alignment_boundaries", "reference_alignment_boundaries", "reference_length", "queryLength", "identities_fraction", "identities_percentage", "matched_fraction", "matched_percentage", "gaps_fraction", "gaps_percentage", "query_length_normalized_2tmscore", "reference_length_normalized_2tmscore", "d0_reference", "d0_query"
    ]
    data_lines = []
    query_length = 0

    alignment_info = {}

    for line in lines:
        # Extract query length from the header
        if "Query (" in line:
            query_length_match = re.search(r"\((\d+) residues\)", line)
            if query_length_match:
                query_length = int(query_length_match.group(1))
        # Extract data lines from the table
        elif re.match(r"^\s*\d+\s+\.\.\.db_structures.*", line):
            data_parts = re.split(r"\s+", line.strip())
            if len(data_parts) >= 10:
                reference_description = "{} {} {}".format(data_parts[1], data_parts[2], data_parts[3])
                data_line = [
                    data_parts[0], reference_description.strip(),
                    data_parts[4], data_parts[5], data_parts[6], data_parts[7],
                    data_parts[8], data_parts[9], data_parts[10], str(query_length)
                ]
                data_lines.append(data_line)
        # Extract alignment details lower in the file
        elif re.match(r"^\d+\.\s*$", line):
            current_no = line.strip().strip('.')
        elif "Identities" in line:
            identities_match = re.search(r"Identities = (\d+/\d+) \((\d+%)\), Matched = (\d+/\d+) \((\d+%)\), Gaps = (\d+/\d+) \((\d+%)\)", line)
            if identities_match:
                if current_no not in alignment_info:
                    alignment_info[current_no] = {}
                alignment_info[current_no].update({
                    "identities_fraction": identities_match.group(1),
                    "identities_percentage": identities_match.group(2),
                    "matched_fraction": identities_match.group(3),
                    "matched_percentage": identities_match.group(4),
                    "gaps_fraction": identities_match.group(5),
                    "gaps_percentage": identities_match.group(6)
                })
        elif "2TM-score" in line:
            twotm_match = re.search(r"2TM-score \(Refn./Query\) = ([\d\.]+) / ([\d\.]+)", line)
            if twotm_match:
                if current_no not in alignment_info:
                    alignment_info[current_no] = {}
                alignment_info[current_no].update({
                    "query_length_normalized_2tm": twotm_match.group(2),
                    "reference_length_normalized_2tm": twotm_match.group(1)
                })
        elif "d0" in line:
            d0_match = re.search(r"d0 \(Refn./Query\) = ([\d\.]+) / ([\d\.]+)", line)
            if d0_match:
                if current_no not in alignment_info:
                    alignment_info[current_no] = {}
                alignment_info[current_no].update({
                    "d0_reference": d0_match.group(1),
                    "d0_query": d0_match.group(2)
                })

    # Add alignment information to data lines
    for data_line in data_lines:
        no = data_line[0]
        if no in alignment_info:
            data_line.extend([
                alignment_info[no].get("identities_fraction", ""),
                alignment_info[no].get("identities_percentage", ""),
                alignment_info[no].get("matched_fraction", ""),
                alignment_info[no].get("matched_percentage", ""),
                alignment_info[no].get("gaps_fraction", ""),
                alignment_info[no].get("gaps_percentage", ""),
                alignment_info[no].get("query_length_normalized_2tm", ""),
                alignment_info[no].get("reference_length_normalized_2tm", ""),
                alignment_info[no].get("d0_reference", ""),
                alignment_info[no].get("d0_query", "")
            ])
        else:
            data_line.extend(["", "", "", "", "", "", "", "", "", ""])

    if not data_lines:
        raise ValueError("Could not find valid data lines in the input file. Please ensure the file contains the expected data format.")

    # Write the output TSV file
    with open(output_file, 'w', newline='') as outfile:
        tsv_writer = csv.writer(outfile, delimiter='\t')
        tsv_writer.writerow(header_line)
        for data_line in data_lines:
            tsv_writer.writerow(data_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a gtalign output file into TSV format.")
    parser.add_argument('--input', required=True, help="Path to the input gtalign output file.")
    parser.add_argument('--output', required=True, help="Path to the output TSV file.")

    args = parser.parse_args()
    parse_gtalign_output(args.input, args.output)
