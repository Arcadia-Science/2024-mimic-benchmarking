import argparse
import os
import subprocess


def decompress_files(input_txt, zip_file, dest_dir):
    with open(input_txt) as f:
        files_to_extract = [line.strip() for line in f if line.strip()]

    os.makedirs(dest_dir, exist_ok=True)

    for file_path in files_to_extract:
        try:
            # Extract each file individually to the specified destination directory
            subprocess.run(
                ["unzip", "-j", zip_file, file_path, "-d", dest_dir],
                check=True,
                stdout=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError as e:
            print(f"Error extracting {file_path}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Decompress specific files from a zip archive based on a file list."
    )
    parser.add_argument(
        "--input_txt",
        type=str,
        required=True,
        help="Path to text file containing file paths to extract from zip archive.",
    )
    parser.add_argument(
        "--zip_file",
        type=str,
        required=True,
        help="Path to zip file containing files to decompress.",
    )
    parser.add_argument(
        "--dest_dir",
        type=str,
        required=True,
        help="Destination directory to store decompressed files.",
    )

    args = parser.parse_args()

    decompress_files(args.input_txt, args.zip_file, args.dest_dir)


if __name__ == "__main__":
    main()
