file_in = snakemake.input[0]
file_out = snakemake.output[0]

with open(file_in, 'r') as f_in, open(file_out, 'w') as f_out:
    lines = f_in.readlines()
    new_lines = []
    for line in lines:
        if line.startswith('>'):
            line_new = line.split('_')
            if len(line_new[0]) < 8:
                new_lines.append(line_new[0] + '_' + line_new[1]+'\n')
            else:   
                new_lines.append(line_new[0]+'\n')
        else:
            new_lines.append(line)
    f_out.writelines(new_lines)



