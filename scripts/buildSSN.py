from foxizyme import ssn, alignment
import argparse

# ## Define input variables
parser = argparse.ArgumentParser()
parser.add_argument('galaxy_txt', help='Input TXT file with galaxy output.')
parser.add_argument('database_fasta', help='Input database fasta')
parser.add_argument('--evalue_threshold', default=1e-20, help='Max allowed evalue.')
parser.add_argument('--cdhit_threshold', default=0.7, help='CD-HIT threshold')
parser.add_argument('--pid_threshold', default=0.5, help='PID SSN threshold to write network file')
parser.add_argument('--fasta_file', default=None, help='Fasta file with sequences to be included explicitly')
parser.add_argument('--target_sequence', default=None, help='Target sequence for residue position extraction.')
parser.add_argument('--target_positions', default=None, help='Comma-separated sequence positions to get variability.')
parser.add_argument('--node_info_file', default='node_info.txt', help='File to write node info')
parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite previous output files')
args=parser.parse_args()

# Store variables
galaxy_txt = args.galaxy_txt
database_fasta = args.database_fasta
evalue_threshold = float(args.evalue_threshold)
cdhit_threshold = float(args.cdhit_threshold)
pid_threshold = float(args.pid_threshold)
node_info_file = args.node_info_file
target_sequence = args.target_sequence
target_positions = args.target_positions
if target_positions != None:
    target_positions = [int(p) for p in target_positions.split(',')]
fasta_file = args.fasta_file
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

if fasta_file != None:
    exp_sequences = alignment.readFastaFile(fasta_file)
    print('Loading %s sequences from fasta file: %s' % (len(exp_sequences), fasta_file))
    for s in exp_sequences:
        cluster_sequences[s] = exp_sequences[s]

print('Building a sequence similarity network for %s sequences' % len(cluster_sequences))

ssn = ssn.sequenceSimilarityNetwork(cluster_sequences, overwrite=overwrite)

sif_file = 'network_'+str('%.2f' % pid_threshold)+'.sif'
print('Writing network sif file %s at a %s PID threshold' % (sif_file, pid_threshold))

ssn.createSifNetworkFile(pid_threshold, overwrite=True, output_file=sif_file)

if target_positions != None and target_sequence == None:
    raise ValueError('You asked for the variability at some target positions, but no target_sequence was given!')

aa_identities = {}
zf = len(str(max(target_positions)))
for t in target_positions:
    label = 'Position_'+str(t).zfill(zf)
    aa_identities[label] = {}
if target_sequence != None:

    if target_positions == None:
        raise ValueError('You asked to get residue variabilities, but no target_positions were given!')

    print('Calculating MSA for %s sequences at msa.fasta' % len(cluster_sequences))
    msa = alignment.mafft.multipleSequenceAlignment(cluster_sequences, stderr=False, stdout=False)
    alignment.writeMsaToFastaFile(msa, 'msa.fasta')

    print('Getting identities of positions %s aligned to %s' % (str(target_positions), target_sequence))
    msa_indexes = alignment.msaIndexesFromSequencePositions(msa, target_sequence, target_positions)
    aas = alignment.getSequencePositionsFromMSAindexes(msa, msa_indexes, return_identity=True)
    for i,label in enumerate(aa_identities):
        for s in aas:
            aa_identities[label][s] = aas[s][i]

    for p in aa_identities:
        ssn.addNodeAttribute(p, aa_identities[p], overwrite=True)

print('Adding sequence length as node_info attribute')
sequence_length = {}
for s in cluster_sequences:
    sequence_length[s] = len(cluster_sequences[s])
ssn.addNodeAttribute('Sequence_length', sequence_length, overwrite=True)

print('Writing node information file at %s' % node_info_file)
ssn.createNodeAttributesFile(node_info_file, overwrite=True)
