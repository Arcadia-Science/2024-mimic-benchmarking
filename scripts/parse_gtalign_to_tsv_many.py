import argparse
import csv
import re
import os
import glob


def parse_gtalign_output(input_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Parse the header and data sections
    header_line = [
        "num", "query", "target", "qlen", "tlen", "alnlen", "qtmscore", "ttmscore", "qcov", "tcov", "qstart", "qend", "tstart", "tend", "rmsd", "fident", "pident", "fmatched", "pmatched", "fgaps", "pgaps", "q2tmscore", "t2tmscore", "td0", "qd0"
    ]
    data_lines = []
    query_length = 0

    alignment_info = {}

    #for line in lines:
    for i, line in enumerate(lines):
        # Extract query length from the header.
        if "Query (" in line:
            query_length_match = re.search(r"\((\d+) residues\)", line)
            if query_length_match:
                query_length = int(query_length_match.group(1))
            if i + 1 < len(lines):
                next_line = lines[i + 1]
                query = str(os.path.basename(next_line))

        # Extract data lines from the table
        # This is a really fragile regex. Might be worth trying to change.
        # It currently matches to a space starting the line (which I think will fail above 99,999),
        # any number of digits, a space, and then a string containing "pdb" (so will only work for pdb files). 
        elif re.match(r"^\s+\d+\s.*pdb", line):
            data_parts = re.split(r"\s+", line.strip())
            if len(data_parts) >= 10:
                target = str(os.path.basename(data_parts[1]))
                qstart, qend = map(int, data_parts[8].split("-"))
                tstart, tend = map(int, data_parts[9].split("-")) 
                qcov = int(data_parts[7]) / query_length
                tcov = int(data_parts[7]) / int(data_parts[10])
                data_line = [
                    data_parts[0], # num
                    query.strip(),
                    target.strip(),
                    str(query_length), # qlen
                    data_parts[10], # tlen
                    data_parts[7], # alnlen
                    data_parts[5], # qtmscore
                    data_parts[4], # ttmscore
                    qcov, tcov, qstart, qend, tstart, tend,
                    data_parts[6], # rmsd
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
                    "fident": identities_match.group(1),
                    "pident": identities_match.group(2).rstrip('%'),
                    "fmatched": identities_match.group(3),
                    "pmatched": identities_match.group(4).rstrip('%'),
                    "fgaps": identities_match.group(5),
                    "pgaps": identities_match.group(6).rstrip('%')
                })
        elif "2TM-score" in line:
            twotm_match = re.search(r"2TM-score \(Refn./Query\) = ([\d\.]+) / ([\d\.]+)", line)
            if twotm_match:
                if current_no not in alignment_info:
                    alignment_info[current_no] = {}
                alignment_info[current_no].update({
                    "q2tm": twotm_match.group(2),
                    "t2tm": twotm_match.group(1)
                })
        elif "d0" in line:
            d0_match = re.search(r"d0 \(Refn./Query\) = ([\d\.]+) / ([\d\.]+)", line)
            if d0_match:
                if current_no not in alignment_info:
                    alignment_info[current_no] = {}
                alignment_info[current_no].update({
                    "td0": d0_match.group(1),
                    "qd0": d0_match.group(2)
                })

    # Add alignment information to data lines
    for data_line in data_lines:
        no = data_line[0]
        if no in alignment_info:
            data_line.extend([
                alignment_info[no].get("fident", ""),
                alignment_info[no].get("pident", ""),
                alignment_info[no].get("fmatched", ""),
                alignment_info[no].get("pmatched", ""),
                alignment_info[no].get("fgaps", ""),
                alignment_info[no].get("pgaps", ""),
                alignment_info[no].get("q2tm", ""),
                alignment_info[no].get("t2tm", ""),
                alignment_info[no].get("td0", ""),
                alignment_info[no].get("qd0", "")
            ])
        else:
            data_line.extend(["", "", "", "", "", "", "", "", "", ""])

    if not data_lines:
        raise ValueError("Could not find valid data lines in the input file. Please ensure the file contains the expected data format.")

    return header_line, data_lines

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse multiple gtalign output files in a directory into a single TSV format.")
    parser.add_argument('--input_dir', required=True, help="Path to the directory containing gtalign output files.")
    parser.add_argument('--output', required=True, help="Path to the combined output TSV file.")

    args = parser.parse_args()

    all_data = []
    header = None

    for file_path in glob.glob(os.path.join(args.input_dir, "*.out")):
        file_header, file_data = parse_gtalign_output(file_path)
        if header is None:
            header = file_header
        all_data.extend(file_data)

    if not all_data:
        raise ValueError("No valid data found in the input directory.")

    # Write the combined output TSV file
    with open(args.output, 'w', newline='') as outfile:
        tsv_writer = csv.writer(outfile, delimiter='\t')
        tsv_writer.writerow(header)
        tsv_writer.writerows(all_data)
