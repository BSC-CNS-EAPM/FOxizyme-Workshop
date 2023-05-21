from foxizyme import ssn, alignment
import argparse

# ## Define input variables
parser = argparse.ArgumentParser()
parser.add_argument('galaxy_txt', help='Input TXT file with galaxy output.')
parser.add_argument('database_fasta', help='Input database fasta')
parser.add_argument('--evalue_threshold', default=1e-20, help='Max allowed evalue.')
parser.add_argument('--cdhit_threshold', default=0.7, help='CD-HIT threshold')
parser.add_argument('--pid_threshold', default=0.5, help='PID SSN threshold to write network file')
parser.add_argument('--calculate_msa', action='store_true', default=False, help='Calculate MSA for SSN sequences?')
parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite previous output files')

args=parser.parse_args()

# Store variables
galaxy_txt = args.galaxy_txt
database_fasta = args.database_fasta
evalue_threshold = float(args.evalue_threshold)
cdhit_threshold = float(args.cdhit_threshold)
pid_threshold = float(args.pid_threshold)
calculate_msa = args.calculate_msa
overwrite = args.overwrite

# Parse galaxy codes
print('Filtering sequences at a maximum e-value of %s' % evalue_threshold)
up_codes = []
leftout = 0
with open(galaxy_txt) as gf:
    for l in gf:
        if 'pid' in l:
            ls = l.split('|')
            if 'gb' in l:
                gid = ls.index('gb')
                code = ls[gid+1]
                e_value = float(ls[-1].split()[4])
                if e_value > evalue_threshold:
                    leftout += 1
                    continue
                up_codes.append(code)
            elif 'pdb' in l:
                pdb = ls.index('pdb')
                code = ls[pdb+1]
                e_value = float(ls[-1].split()[4])
                if e_value > evalue_threshold:
                    leftout += 1
                    continue
                up_codes.append(code)

print('\t%s sequences were selected.' % len(up_codes))
print('\t%s sequences were left out.' % leftout)

# Parse database fasta
print("Reading %s sequences from the database: %s" % (len(up_codes), database_fasta))
database_sequences = alignment.readFastaFile(database_fasta)
sequences = {}
for s in database_sequences:
    ss = s.split('|')
    if 'gb' in s:
        gid = ss.index('gb')
        code = ss[gid+1]
        # if code.endswith('.1'):
        #     code = code[:-2]

        sequences[code] = database_sequences[s]
    elif 'pdb' in s:
        pdb = ss.index('pdb')
        code = ss[pdb+1]
        sequences[code] = database_sequences[s]

target_sequences = {}
for code in up_codes:
    try:
        target_sequences[code] = sequences[code]
    except:
        # print('Sequence %s was not found in the database' % code)
        continue

# Clustering
print('Clustering sequences with CD-HIT at %s PID threshold' % cdhit_threshold)
clusters = alignment.cdhit.clusterSequences(target_sequences, pid_threshold=cdhit_threshold)

print('Found %s clusters.' % len(clusters))

cluster_sequences = {}
for c in clusters:
    cluster_sequences[c] = target_sequences[c]

print('Building a sequence similarity network for %s sequences' % len(cluster_sequences))

ssn = ssn.sequenceSimilarityNetwork(cluster_sequences, overwrite=overwrite)

sif_file = 'network_'+str('%.2f' % pid_threshold)+'.sif'
print('Writing network sif file %s at a %s PID threshold' % (sif_file, pid_threshold))

ssn.createSifNetworkFile(pid_threshold, overwrite=True, output_file=sif_file)

if calculate_msa:
    print('Calculating MSA for %s sequences at msa.fasta' % len(cluster_sequences))
    msa = alignment.mafft.multipleSequenceAlignment(cluster_sequences, stderr=False, stdout=False)
    alignment.writeMsaToFastaFile(msa, 'msa.fasta')
