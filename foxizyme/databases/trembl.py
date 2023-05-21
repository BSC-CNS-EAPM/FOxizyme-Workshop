from .. import alignment
import os

def BLAST(job_folder, target_sequences, num_descriptions=500,
         trembl_path='/gpfs/projects/bsc72/uniprot/trembl/trembl'):
    """
    Do a BLAST search in the TREMBL database at a BSC cluster.

    Parameters
    ==========
    target_sequences : dict
        The sequences for which to run BLAST against TREMBL.
    num_descriptions : int
        The maximum number of outputs to get in the BLAST output.
    """
    if not isinstance(target_sequences, dict):
        raise ValueError('target_sequences must be a dictionary!')

    # Creat job folders
    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Iterate each given sequence
    jobs = []
    for sequence in target_sequences:

        # Create sequence folder
        sdir = job_folder+'/'+sequence
        if not os.path.exists(sdir):
            os.mkdir(sdir)

        # Save fasta files in inputs folder
        alignment.writeFastaFile({sequence:target_sequences[sequence]},
                                 sdir+'/'+sequence+'.fasta')

        # Create command line for running PSI BLAST
        command = 'cd '+sdir+'\n'
        command += 'blastp '
        command += '-query '+sequence+'.fasta '
        command += '-num_descriptions '+str(num_descriptions)+' '
        command += '-db '+trembl_path+' '
        command += '-out '+sequence+'.out '
        command += '\n'
        command += 'cd ../../\n'

        jobs.append(command)

    return jobs

def PSIBlast(job_folder, target_sequences, num_descriptions=5000, num_iterations=10, evalue=10, max_target_seqs=None, output_pssm=False):
    """
    Do a PSIBLAST search in the TREMBL database at a BSC cluster.

    Parameters
    ==========
    target_sequences : dict

    """
    if not isinstance(target_sequences, dict):
        raise ValueError('target_sequences must be a dictionary!')

    trembl_path = '/gpfs/projects/bsc72/uniprot/trembl/trembl'

    # Creat job folders
    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Iterate each given sequence
    jobs = []
    for sequence in target_sequences:

        # Create sequence folder
        sdir = job_folder+'/'+sequence
        if not os.path.exists(sdir):
            os.mkdir(sdir)

        # Save fasta files in inputs folder
        alignment.writeFastaFile({sequence:target_sequences[sequence]},
                                 sdir+'/'+sequence+'.fasta')

        if num_descriptions != None and max_target_seqs != None:
            raise ValueError('Using num_descriptions and max_target_seqs is not compatible. Choose one of the two.')

        # Create command line for running PSI BLAST
        command = 'cd '+sdir+'\n'
        command += 'psiblast '
        command += '-query '+sequence+'.fasta '
        if num_descriptions != None:
            command += '-num_descriptions '+str(num_descriptions)+' '
        if max_target_seqs != None:
            command += '-max_target_seqs '+str(max_target_seqs)+' '
        command += '-db '+trembl_path+' '
        command += '-out '+sequence+'.out '
        if output_pssm:
            command += '-out_ascii_pssm pssm.smp -save_each_pssm '
        command += '-num_iterations '+str(num_iterations)+' '
        command += '-evalue '+str(evalue)
        command += '\n'
        command += 'cd ../../\n'

        jobs.append(command)

    return jobs

def readBLASTResults(job_folder, only_codes=True):
    """
    Read the resuts from the PSIBlast() function over the TREMBL database.

    Parameters
    ==========
    job_folder : str
        Path to the job folder where the output of the PSIBlast() calculation is contained.
    only_codes : bool
        Split the sequence headers to only get the Uniprot codes.

    Returns
    =======
    sequences : dict
        The results from the PSIBLAST calculation.
    """

    sequences = {}
    for model in os.listdir(job_folder):

        model_out = job_folder+'/'+model+'/'+model+'.out'
        if not os.path.exists(model_out):
            print('BLAST output for model %s was not found!')
            continue
        else:
            blast_results = alignment.blast._parseBlastpOutput(model_out)
            sequences[model] = {}
            for s in blast_results:
                if only_codes:
                    sequences[model][s.split()[0]] = blast_results[s]['e-value']
                else:
                        sequences[model].append(s)

    return sequences

def readPSIBLASTResults(job_folder, as_one_bundle=False, only_codes=False):
    """
    Read the resuts from the PSIBlast() function over the TREMBL database.

    Parameters
    ==========
    job_folder : str
        Path to the job folder where the output of the PSIBlast() calculation is contained.
    as_one_bundle : bool
        Get the sequences as one unique list, i.e., without separation for iterations.
    only_codes : bool
        Split the sequence headers to only get the Uniprot codes.

    Returns
    =======
    sequences : dict
        The results from the PSIBLAST calculation.
    """

    sequences = {}
    for model in os.listdir(job_folder):

        model_out = job_folder+'/'+model+'/'+model+'.out'
        if not os.path.exists(model_out):
            print('PSIBLAST output for model %s was not found!' % model)
            continue
        else:
            psiblast_results = alignment.blast._parsePSIBlastOutput(model_out)
            if as_one_bundle:
                sequences[model] = {}
                for c in psiblast_results:
                    for s in psiblast_results[c]:
                        if only_codes:
                            sequences[model][s.split()[0]] = psiblast_results[c][s]['e-value']
                        else:
                            sequences[model][s] = psiblast_results[c][s]['e-value']
            else:
                sequences.setdefault(model, {})
                #sequences[model] = {}

                for c in psiblast_results:
                    sequences[model][c] = {}
                    for s in psiblast_results[c]:
                        if only_codes:
                            sequences[model][c][s.split()[0]] = psiblast_results[c][s]['e-value']
                        else:
                            sequences[model][c][s] = psiblast_results[c][s]['e-value']

        if as_one_bundle:
            sequences[model] = list(set(sequences[model]))

    return sequences
