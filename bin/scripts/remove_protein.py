# Remove protein identifiers from FASTA file and only keeps plasmid name
# Input: FASTA file with protein sequence
# Output: FASTA file with only the header lines (without protein sequence)

# Read in fasta file
file_in = snakemake.input[0]

# Write output to a new fasta file
file_out = snakemake.output[0]

print(f"Reading from {file_in}...")
with open(file_in, 'r') as f_in, open(file_out, 'w') as f_out:
    # Read in all lines from the input file
    lines = f_in.readlines()

    print(f"Read {len(lines)} lines from input file.")

    # Initialize a new list to store the modified lines
    new_lines = []

    # Iterate over each line in the file
    for i, line in enumerate(lines):
        print(f"Processing line {i+1}...")

        # Check if the line starts with '>'
        if line.startswith('>'):
            # Split the line on '_'
            line_new = line.split('_')

            # If the first part of the split line is less than 8 characters long
            if len(line_new[0]) < 8:
                # Add the first and second parts of the split line back together,
                # with an underscore in between, and a newline at the end
                new_line = line_new[0] + '_' + line_new[1] + '\n'
                print(f"Line {i+1}: Writing '{new_line}' to output file.")
                new_lines.append(new_line)
            # Otherwise, just add the first part of the split line, with a newline at the end
            else:
                new_line = line_new[0] + '\n'
                print(f"Line {i+1}: Writing '{new_line}' to output file.")
                new_lines.append(new_line)

        # If the line doesn't start with '>', just add it to the new list as is
        else:
            print(f"Line {i+1}: Writing '{line}' to output file.")
            new_lines.append(line)

    # Write all the new lines to the output file
    print(f"Writing {len(new_lines)} lines to output file.")
    f_out.writelines(new_lines)