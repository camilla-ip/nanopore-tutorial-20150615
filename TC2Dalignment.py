#!/usr/bin/env python

# =========================================================================== #
# TC2Dalignment.py                                                            #
#                                                                             #
# Extract the alignment of the template, complement and 2D read.              #
#                                                                             #
# Camilla Ip                                                                  #
# camilla.ip@well.ox.ac.uk                                                    #
# 2015-06-15                                                                  #
# =========================================================================== #

import h5py, itertools, numpy as np, os, sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# =========================================================================== #
# Command-line.
# =========================================================================== #

if len(sys.argv) < 4:
    print 'Usage: TC2Dalignment.py fast5file outdir debug0or1'
    print '       Extract the alignment of the template to the complement, and the corresponding 2D read.'
    print '       Extract a single fasta file containing the template, complement and 2D nucleotide sequence alignment,'
    print '       with dashes for gaps. If there is no 2D or complement read, only the template is printed.'
    sys.exit(1)
_fast5path = os.path.expandvars(sys.argv[1])
_outdir = os.path.expandvars(sys.argv[2])
_debug = int(sys.argv[3])

# =========================================================================== #
# Constants.
# =========================================================================== #

_linewidth = 100
_ErrorInvalidData = 18

# =========================================================================== #
# Global variables.
# =========================================================================== #

_fast5file = os.path.basename(_fast5path)
_fast5stem = ''.join(_fast5file.split('.')[0:-1])

_R = {}		# Dictionary of raw data from the basecalled fast5 file.
		# keys: [temp|comp|twod]_[events|seq]
                # If the data is not present in the file or could not be extracted, the key exists but the value is None.

_anchor = {}	# Dictionary of anchor points for the temp, comp and twod sequences.
		# _anchor['temp'] = [ [temp_idx, temp_base, model_kmer], ... ]
		# _anchor['comp'] = [ [comp_idx, comp_base, model_kmer], ... ]
		# _anchor['twod'] = [ [temp_idx, comp_idx, twod_base, twod_kmer], ... ]
		# where
                #  the indices are the 0-based indicies into the respective events or alignment table in the fast5 file,
                #  negative indices correspond to fillers that were added to make the sequences match up exactly with the
                #    fastq sequences from the ONT basecaller,
                #  upper case bases were transcribed exactly from the Events table,
                #  lower-case were inferred from either the fastq sequence computed by ONT basecaller.

# =========================================================================== #
# Functions
# =========================================================================== #

def numoffset(s1, s2):
    if s1 == s2:
        return 0
    length = len(s1)
    for i in range(1, length):
        if s1[i:] == s2[:-i]:
            return i
    return length

def print_fastaseq(recordid, seq, filename, linewidth, outfp):
    'Print the fasta record in linewidth-chunks. seq must be a string.'
    seqformatted = "\n".join([seq[i:i+int(linewidth)] for i in range(0, len(seq), int(linewidth))])
    outfp.write('>{0} len={1} file={2}\n{3}\n'.format(recordid, len(seq), filename, seqformatted))

def Extract_RawDataFromFast5(fast5path):
    'Extract fast5 data needed to reconstruct the alignment, return a dictionary of numpy arrays.'

    D = {}
    keys = [ 'temp_events', 'comp_events', 'twod_aln', 'temp_seq', 'comp_seq', 'twod_seq' ]
    for key in keys:
        D[key] = None
    failed = False
    failed_msg = []
    try:
        hdf = h5py.File(fast5path, 'r')
    except:
        failed_msg.append('Failed to access fast5 file')
        failed = True
    try:
        D['channel_id'] = hdf['UniqueGlobalKey/channel_id'].attrs['channel_number']
        D['read_id'] = [x for x in hdf['Analyses/EventDetection_000/Reads'] if x.startswith('Read_')][0].split('_')[1]
        D['file_id'] = os.path.basename(_fast5path).split('_')[-2].replace('file', '')
        D['run_number'] = '_'.join(_fast5file.split('_')[-5:-3])
    except:
        failed = True
    try:
        D['temp_events'] = hdf['Analyses/Basecall_2D_000/BaseCalled_template/Events'][()]
    except:
        failed_msg.append('Failed to retrieve template Events table')
        failed = True
    try:
        D['comp_events'] = hdf['Analyses/Basecall_2D_000/BaseCalled_complement/Events'][()]
    except:
        failed_msg.append('Failed to retrieve complement Events table')
        failed = True
    try:
        D['twod_aln'] = hdf['Analyses/Basecall_2D_000/BaseCalled_2D/Alignment'][()]
    except:
        failed_msg.append('Failed to retrieve 2D alignment')
        failed = True
    try:
        D['temp_seq'] = hdf['Analyses/Basecall_2D_000/BaseCalled_template/Fastq'][()].split('\n')[1]
    except:
        failed_msg.append('Failed to retrieve template fastq string')
        failed = True
    try:
        D['comp_seq'] = hdf['Analyses/Basecall_2D_000/BaseCalled_complement/Fastq'][()].split('\n')[1]
    except:
        failed_msg.append('Failed to retrieve complement fastq string')
        failed = True
    try:
        D['twod_seq'] = hdf['Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'][()].split('\n')[1]
    except:
        failed_msg.append('Failed to retrieve 2D fastq string')
        failed = True
    try:
        hdf.close()
    except:
        pass
    if failed:
        #if _debug:
        #    logpath = os.path.join(_outdir, _fast5stem+'_aln.log')
        #    with open(logpath, 'w') as log_fp:
        #        log_fp.write('{0}\n'.format(
        #            '\n'.join(['Error ({0}): {1}'.format(_fast5file, x) for x in failed_msg])))
        sys.exit(_ErrorInvalidData)
    return D

def Compute_1DAnchor(temp_events, datatype):
    '''
    Collate the 1D anchor points and return as a list of lists.
    How do the template events relate to the fastq sequence?
    The template sequence in the fastq file is constructed from the template Events by
    prefixing with 'GC', then taking the middle base of the 'model_state' for each
    event with move == 1, ignoring all events with move == 0, and inserting (move-1)
    gap characters (which are filled by a process I don't understand in the fastq seq)
    before all events with move > 1.
    '''

    global _anchor

    # If no data of this type in the fast5 file, return.
    if _R[datatype+'_events'] is None:
        return

    # Dump the events table to a file.
    if _debug:
        out_path = os.path.join(_outdir, _fast5stem+'_'+datatype+'_events.txt')
        with open(out_path, 'w') as out_fp:
            H = ['mean', 'start', 'stdv', 'length', 'model_state',
                 'model_level', 'move', 'p_model_state', 'mp_state', 'p_mp_state',
                 'p_A', 'p_C', 'p_G', 'p_T']
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in H])))
            for e in _R[datatype+'_events']:
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in e])))

    # Collate the anchor points and sequence for the 1D sequences.
    _anchor[datatype] = []
    if datatype == 'temp':
        _anchor[datatype].append( [-2, 'g', '-----'] )
        _anchor[datatype].append( [-1, 'c', '-----'] )
    elif datatype == 'comp':
        _anchor[datatype].append( [-2, 't', '-----'] )
        _anchor[datatype].append( [-1, 'c', '-----'] )
    i = 0	# Counter for Event table records.
    j = 2	# Counter for bases in fastq sequence from ONT.
    for e in _R[datatype+'_events']:
        # Skip any Event record where move == 0
        if e[6] == 0 and i > 0:
            i += 1
            continue
        # When move > 1, insert the corresponding bases from the fastq seq called by ONT.
        if e[6] > 1:
            for k in range(1, e[6]):
                _anchor[datatype].append( [float('{0}.{1}'.format(i-1, k)), _R[datatype+'_seq'][j].lower(), '-----'] )
                j += 1
        _anchor[datatype].append( [i, e[4][2], e[4]] )
        i += 1
        j += 1
    if datatype == 'temp':
        _anchor[datatype].append( [i, 'c', '-----'] )
        _anchor[datatype].append( [i+1, 'g', '-----'] )
    elif datatype == 'comp':
        _anchor[datatype].append( [i, 'a', '-----'] )
        _anchor[datatype].append( [i+1, 't', '-----'] )

    # Save the anchor table to a .txt file for checking later.
    if _debug:
        anchor_path = os.path.join(_outdir, _fast5stem+'_'+datatype+'_anchor.txt')
        with open(anchor_path, 'w') as out_fp:
            H = [datatype+'_idx', datatype+'_base']
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in H])))
            for e in _anchor[datatype]:
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in e])))

    # Print error message if inferred sequence is not the same as the ONT sequence.
    inferred_seq = ''.join([x[1] for x in _anchor[datatype]]).upper()
    ont_seq = _R[datatype+'_seq'].upper()
    if inferred_seq != ont_seq:
        positions_that_differ = []
        for i in range(0, min(len(inferred_seq), len(ont_seq))):
            if inferred_seq[i] != ont_seq[i]:
                positions_that_differ.append(i+1)
        if len(inferred_seq) != len(ont_seq):
            for i in range(min(len(inferred_seq), len(ont_seq)), max(len(inferred_seq), len(ont_seq))):
                positions_that_differ.append(i+1)
        sys.stderr.write('Error: Inferred {0} sequence disagrees with ONT sequence at positions {1}\n'.format(
            datatype, ','.join([str(x) for x in positions_that_differ])))

    # Save the predicted and actual sequences to a .fasta file for checking later.
    if _debug:
        fasta_path = os.path.join(_outdir, _fast5stem+'_'+datatype+'_checks.fasta')
        with open(fasta_path, 'w') as out_fp:
            print_fastaseq(datatype+'_pred', ''.join([x[1] for x in _anchor[datatype]]), _fast5file, _linewidth, out_fp)
            print_fastaseq(datatype+'_seq', ''.join(_R[datatype+'_seq']), _fast5file, _linewidth, out_fp)

def Compute_2DAnchor(twod_aln):
    'Collate the 2D anchor information, saving to debug files.'

    # If no data of this type in the fast5 file, return.
    if _R['twod_aln'] is None:
        return

    if _debug:
        twodaln_path = os.path.join(_outdir, _fast5stem+'_'+'twod_alignment.txt')
        with open(twodaln_path, 'w') as out_fp:
            H = ['template', 'complement', 'kmer']
            out_fp.write('{0}\n'.format('\t'.join(H)))
            for e in _R['twod_aln']:
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in e])))

    # Collate the anchor points for the aligment of the template and the complement to create the 2D sequence.
    kmer_prev = 'N'
    _anchor['twod'] = []
    _anchor['twod'].append( [_R['twod_aln'][0][0]-2, _R['twod_aln'][0][1]+2, 'g', '-----'] )
    _anchor['twod'].append( [_R['twod_aln'][0][0]-1, _R['twod_aln'][0][1]+1, 'c', '-----'] )
    i = 0	# Counter for 2D alignment table records.
    j = 2	# Counter for bases in 2D fastq sequence from ONT.
    for e in _R['twod_aln']:
        template, complement, kmer = e
        steplen = numoffset(kmer_prev, kmer)
        if steplen == 0:
            pass
        else:
            if steplen > 1:	# Fill missing 2D bases in alignment with those from the 2D seq from ONT.
                for k in range(1, steplen):
                    _anchor['twod'].append( [-1, -1, _R['twod_seq'][j+k-1].lower(), '-----'] )
                    j += 1
            _anchor['twod'].append( [template, complement, kmer[2], kmer] )
            j += 1
        i += 1
        kmer_prev = kmer
    _anchor['twod'].append( [_anchor['twod'][-1][0]+1, _anchor['twod'][-1][1]-1, 't', '-----'] )
    _anchor['twod'].append( [_anchor['twod'][-1][0]+1, _anchor['twod'][-1][1]-1, 'a', '-----'] )
    # Debug: print '\n'.join(['\t'.join([str(y) for y in x]) for x in _anchor['twod']])

    # Save the anchor table to a .txt file for checking later.
    if _debug:
        anchor_path = os.path.join(_outdir, _fast5stem+'_'+'twod_anchor.txt')
        with open(anchor_path, 'w') as out_fp:
            H = ['temp_idx', 'comp_idx', 'twod_base', 'twod_kmer']
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in H])))
            i = 0
            for e in _anchor['twod']:
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in e])))
                i += 1

    # Print error message if inferred sequence is not the same as the ONT sequence.
    inferred_seq = ''.join([x[2] for x in _anchor['twod']]).upper()
    ont_seq = _R['twod_seq'].upper()
    if inferred_seq != ont_seq:
        positions_that_differ = []
        for i in range(0, min(len(inferred_seq), len(ont_seq))):
            if inferred_seq[i] != ont_seq[i]:
                positions_that_differ.append(i+1)
        if len(inferred_seq) != len(ont_seq):
            for i in range(min(len(inferred_seq), len(ont_seq)), max(len(inferred_seq), len(ont_seq))):
                positions_that_differ.append(i+1)
        sys.stderr.write('Error: Inferred {0} sequence disagrees with ONT sequence at positions {1}\n'.format(
            'twod', ','.join([str(x) for x in positions_that_differ])))

    # Save the predicted and actual sequences to a .fasta file for checking later.
    if _debug:
        fasta_path = os.path.join(_outdir, _fast5stem+'_'+'twod_checks.fasta')
        with open(fasta_path, 'w') as out_fp:
            print_fastaseq('twod_pred', ''.join([x[2] for x in _anchor['twod']]), _fast5file, _linewidth, out_fp)
            print_fastaseq('twod_seq', ''.join(_R['twod_seq']), _fast5file, _linewidth, out_fp)

def Compute_Alignment():
    'Use the _anchor tables to reconstruct the alignment of the template, complement and 2D reads.'

    global _anchor

    if _anchor['temp'] is None or _anchor['comp'] is None or _anchor['twod'] is None:
        return

    # Create dictionaries from the temp and comp anchor data so I can look up values by idx.
    anchorD = {}
    anchorD['temp'] = dict([[x[0], x[1:]] for x in _anchor['temp']])
    anchorD['comp'] = dict([[x[0], x[1:]] for x in _anchor['comp']])
    anchorD['twod'] = dict([[x[0], x[1:]] for x in _anchor['twod']])

    # Collate the alignment of the three sequences, based on the anchor points.
    _anchor['align'] = []		# elt = [temp_idx, temp_base, comp_idx, comp_base, twod_base]
    prev_e = [-3, sys.maxint, 'N', '-----']	# temp_idx, comp_idx, twod_base, twod_kmer
    i = 0
    numanchors = len(_anchor['twod'])
    for e in _anchor['twod']:
        temp_idx, comp_idx, twod_base, twod_kmer = e
        # If the temp_idx has jumped, fill in the missing ones. 
        if (float(e[0]) - float(prev_e[0])) > 1:
            keyL = [x for x in set([float(x[0]) for x in _anchor['temp'] \
                    if float(x[0]) > float(prev_e[0]) and float(x[0]) < float(e[0])] + range(prev_e[0]+1, e[0]))]
            keyL.sort(key=float)
            for idx in keyL:
                if float(int(idx)) == idx:
                    idx = int(idx)
                L = [idx, anchorD['temp'][idx][0] if anchorD['temp'].has_key(idx) else '-', -1, '-', '-']
                _anchor['align'].append(L)
        # If the comp_idx has jumped, fill in the missing ones.
        if (float(prev_e[1]) - float(e[1])) > 1 and prev_e[1] != sys.maxint and e[1] != -1:
            keyL = [x for x in set([float(x[0]) for x in _anchor['comp'] \
                    if float(x[0]) > float(e[1]) and float(x[0]) < float(prev_e[1])] + range(e[1]+1, prev_e[1]))]
            keyL.sort(key=float)
            keyL.reverse()
            for idx in keyL:
                if float(int(idx)) == idx:
                    idx = int(idx)
                L = [-1, '-', idx, anchorD['comp'][idx][0] if anchorD['comp'].has_key(idx) else '-', '-']
                _anchor['align'].append(L)
        # Print the current record.
        temp_base = anchorD['temp'][temp_idx][0] if anchorD['temp'].has_key(temp_idx) and i > 1 and temp_idx >= 0 else '-'
        comp_base = anchorD['comp'][comp_idx][0] if anchorD['comp'].has_key(comp_idx) and i > 1 and comp_idx >= 0 else '-'
        L = [temp_idx, temp_base, comp_idx, comp_base, twod_base]
        _anchor['align'].append(L)
        # Bookkeeping
        if e[0] > prev_e[0]:
            prev_e[0] = e[0]
        if e[1] < prev_e[1] and (e[1] != -1 or (numanchors - i) <= 2):
            prev_e[1] = e[1]
        prev_e[2] = e[2]
        prev_e[3] = e[3]
        i += 1

    # Save the alignment of the three files to a log file.
    out_path = os.path.join(_outdir, _fast5stem+'_'+'all3_anchor.txt')
    with open(out_path, 'w') as out_fp:
            H = ['temp_idx', 'temp_base', 'comp_idx', 'comp_base', 'twod_base']
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in H])))
            for e in _anchor['align']:
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in e])))

    # Save the alignment of the three files to a fasta file.
    fasta_path = os.path.join(_outdir, _fast5stem+'_'+'all3_aln.fasta')
    with open(fasta_path, 'w') as out_fp:

        recordid = '{run_number}_channel_{channel_id}_read_{read_id}_temp'.format(
            run_number=_R['run_number'], channel_id=_R['channel_id'], read_id=_R['read_id'])
        seq = ''.join([x[1] for x in _anchor['align']])
        print_fastaseq(recordid, seq, _fast5file, _linewidth, out_fp)

        recordid = '{run_number}_channel_{channel_id}_read_{read_id}_compc'.format(
            run_number=_R['run_number'], channel_id=_R['channel_id'], read_id=_R['read_id'])
        seq = seq = str(Seq(''.join([x[3] for x in _anchor['align']]), IUPAC.unambiguous_dna).complement())
        print_fastaseq(recordid, seq, _fast5file, _linewidth, out_fp)

        recordid = '{run_number}_channel_{channel_id}_read_{read_id}_2D'.format(
            run_number=_R['run_number'], channel_id=_R['channel_id'], read_id=_R['read_id'])
        seq = ''.join([x[4] for x in _anchor['align']])
        print_fastaseq(recordid, seq, _fast5file, _linewidth, out_fp)

        recordid = '{run_number}_channel_{channel_id}_read_{read_id}_2Danchored'.format(
            run_number=_R['run_number'], channel_id=_R['channel_id'], read_id=_R['read_id'])
        seq = ''.join([(x[4] if x[0] > -1 and x[2] > -1 else ('-' if x[4] == '-' else 'N')) for x in _anchor['align']])
        print_fastaseq(recordid, seq, _fast5file, _linewidth, out_fp)

# =========================================================================== #
# Main
# =========================================================================== #

if __name__ == '__main__':

    _R = Extract_RawDataFromFast5(_fast5path)
    Compute_1DAnchor(_R['temp_events'], 'temp')
    Compute_1DAnchor(_R['comp_events'], 'comp')
    Compute_2DAnchor(_R['twod_aln'])
    Compute_Alignment()
    sys.exit(0)

# =========================================================================== #
