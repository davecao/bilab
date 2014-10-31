# -*- coding: utf-8 -*-

from __future__ import print_function

import sys
import os.path
import time
import datetime
import pprint
import tempfile


from subprocess import Popen, PIPE, check_call
from threading import Thread

try:
    from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty

from bilab.sequence import PSSM
from bilab.utilities import AsynchronousFileReader

__all__ = ['NCBIWWW', 'NCBIStandalone']


def generate_temp_query(id, sequence):
    tmpdir = tempfile.gettempdir()
    temp = tempfile.NamedTemporaryFile(suffix='.fasta',
                                       prefix='deltablast_',
                                       delete=False,
                                       dir=tmpdir)
    try:
        temp.writelines('>' + id + '\n' + sequence)
    finally:
        temp.close()
    return temp.name

def get_num_physical_cpus():
    # use package hwloc
    try:
        import hwloc
        topology = hwloc.Topology()
        topology.load()
        return topology.get_nbojjs_by_type(hwloc.OBJ_CORE)
    except (ImportError, NotImplementError):
        pass

def get_num_logical_cpus():
    # ver 2.6+
    # get logical cpus
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementError):
        pass

def run(cmd, out=PIPE, err=PIPE, verbose=True):
    """ Run command in subprocess """
    if verbose:
        print("Command:\n{0} \n Waiting ... "
                .format(" ".join(str(i) for i in cmd)), end="")
    # launch the command
    proc = Popen(cmd, shell=False, stdout=PIPE, stderr=PIPE)
    #ret = check_call(cmd, shell=False, stdout=PIPE, stderr=PIPE)

    # Launch the asynchronous readers of the process' stdout and stderr
    stdout_queue = Queue()
    stdout_reader = AsynchronousFileReader(proc.stdout, stdout_queue)
    stdout_reader.start()
    stderr_queue = Queue()
    stderr_reader = AsynchronousFileReader(proc.stderr, stderr_queue)
    stderr_reader.start()
    output = ''
     # Check the queues if we received some output (until there is nothing more to get).
    while not stdout_reader.eof() or not stderr_reader.eof():
        # Show what we received from standard output.
        while not stdout_queue.empty():
            line = stdout_queue.get()
            print("stdout:{0}".format(repr(line)))
            output += line.strip()

        # Show what we received from standard error.
        while not stderr_queue.empty():
            line = stderr_queue.get()
            print("stderr: {0}".format(repr(line)))

        # Sleep a bit before asking the readers again.
        time.sleep(.1)

    # Let's be tidy and join the threads we've started.
    stdout_reader.join()
    stderr_reader.join()

    # Close subprocess' file descriptors.
    proc.stdout.close()
    proc.stderr.close()
    if verbose:
        print("Finished.")
    return output

def openURL(url, timeout=5, **kwargs):
    """ Using urlib/urllib2 """
    try:
        from urllib2 import urlopen, URLError, Request
    except ImportError:
        from urllib.request import urlopen, Request
        from urllib.error import URLError
    if kwargs:
        request = Request(url, **kwargs)
    else:
        request = str(url)
    try:
        return urlopen(request, timeout=int(timeout))
    except URLError:
        raise IOError('{0] could not be opened for reading, invalid URL or '
                    'no internet connection'.format(repr(request)))

def clear(last_msg):
    """ Clear console """
    sys.stderr.write('\r' + ' ' * (len(last_msg)) + '\r')

def write(msg):
    """ Write msg to console """
    sys.stderr.write(msg)
    sys.stderr.flush()

def sleep_counter(seconds, msg=''):
    """
        Sleep for seconds while updating screen message every second.
        Message will start with ``'Waiting for Xs '`` followed by *msg*.
    """
    msg = str(msg)
    for second in range(int(seconds), 0, -1):
        info = 'Waiting for {0}s {1}'.format(second, msg)
        write(info)
        time.sleep(1)
        clear(info)

def NCBIWWW(sequence,
            filename=None,
            program = 'blastp',
            blast_programs='deltaBlast',
            db = 'pdb', **kwargs):

#                   db = 'pdb',
#                   evalue = 10e-10,
#                   hitsize = 250,
#                   Threshold=0.005,
#                   Pseudocount=0, **kwargs):
    """
    Interface for standalone version, only support deltaBlast now.

    """
    if sequence == 'test':
        sequence=('MIEIKDKQLTGLRFIDLFAGLGGFRLALESCGAECVYSNEWDKYAQEVYEMNFGEKPEGD'
            'ITQVNEKTIPDHDILCAGFPCQAFSISGKQKGFEDSRGTLFFDIARIVREKKPKVVFMENVKNFAS'
            'HDNGNTLEVVKNTMNELDYSFHAKVLNALDYGIPQKRERIYMICFRNDLNIQNFQFPKPFELNTFV'
            'KDLLLPDSEVEHLVIDRKDLVMTNQEIEQTTPKTVRLGIVGKGGQGERIYSTRGIAITLSAYGGGI'
            'FAKTGGYLVNGKTRKLHPRECARVMGYPDSYKVHPSTSQAYKQFGNSVVINVLQYIAYNIGSSLNF'
            'KPY')
    try:
        sequence = ''.join(sequence.split())
        _ = sequence.isalpha()
    except AttributeError:
        raise ValueError('sequence must be a string.')
    else:
        if not _:
            raise ValueError('not a valid sequence')
    # set query header
    headers = {'User-agent':'bilab'}
    # basic query
    query = [
             ('DATABASE', db),
             ('ENTREZ_QUERY','(none)'),
             ('PROGRAM', program),
             ('BLAST_PROGRAMS', blast_programs)
            ]
    # set expect
    expect = float(kwargs.pop('expect', 10e-10))
    if expect <= 0:
        raise ValueError('expect must be a positive number')
    query.append(('EXPECT', expect))
    # set hitlist size
    hitlist_size = int(kwargs.pop('hitlist_size', 250))
    if hitlist_size <= 0:
        raise ValueError('hitlist_size should be a positive number')
    query.append(('HITLIST_SIZE', hitlist_size))
    query.append(('QUERY', sequence))
    query.append(('CMD', 'Put'))

    sleep = int(kwargs.pop('sleep', 2))
    timeout = float(kwargs.pop('timeout', 120))
    if kwargs:
        print("Keyword arguments {0} are not used"
            .format(", ".join([repr(key) for key in kwargs])))

    try:
        import urllib.parse
        urlencode = lambda data: bytes(urllib.parse.urlencode(data), 'utf-8')
    except ImportError:
        from urllib import urlencode

    url = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi'
    data = urlencode(query)

    print('Blast searching NCBI PDB database for "{0} ... "'
                .format(sequence[:5]))

    handle = openURL(url, data=data, headers=headers)
    html = handle.read()

    index = html.find(b'RID =')
    if index == -1:
        raise Exception('NCBI did not return expected response')
    else:
        last = html.find(b'\n', index)
        rid = html[index + len('RID ='):last].strip()
        #print("Query data: {0}".format(data))
        print("Get request id: {0}".format(rid))

    index = html.find(b'RTOE =')
    if index == -1:
        rtoe = None
    else:
        last = html.find(b'\n', index)
        rtoe = int(html[index + len('RTOE ='):last].strip())

    query = [('ALIGNMENTS', 500), ('DESCRIPTIONS', 500),
             ('FORMAT_TYPE', 'XML'), ('RID', rid), ('CMD', 'Get')]

    data = urlencode(query)

    start_time = time.time()

    while True:
        sleep_counter(sleep, msg = 'to reconnect NCBI for search results...')
        info = 'Connecting NCBI for search results...'
        write(info)
        handle = openURL(url, data=data, headers=headers)
        results = handle.read()
        index = results.find(b'Status=')
        clear(info)
        if index < 0:
            break
        last = results.index(b'\n', index)
        status = results[index+len('Status='):last].strip()
        if status.upper() == 'READY':
            break
        sleep = int(sleep * 1.5)
        elapse_time = time.time() - start_time
        if elapse_time > timeout:
            print('Blast search time out.')
            return None
    print('Blast search completed in %.1fs.', elapse_time)
    try:
        ext_xml = filename.lower().endswith('.xml')
    except AttributeError:
        pass
    else:
        if not ext_xml:
            filename += '.xml'
        out = open(filename, 'w')
        out.write(results)
        out.close()
        print('Results are saved as {0}.'.format(repr(filename)))

def NCBIStandalone(sequence,
                   queryid=None,
                   filename = None,
                   program = 'deltablast',
                   deletetmp = True,
                   verbose = True,
                   **kwargs):
    """
     1. Download cdd_delta.tar.gz data from ftp://ftp.ncbi.nlm.nih.gov/blast/db

     2. tar zxvf cdd_delta.tar.gz into /opt/work/db/cdd_delta
     
     3. Set user's ncbirc

    Example of .ncbirc::

        [NCBI]
        Data=/opt/apps/bioMatrices

        ; Start the section for BLAST configuration
        [BLAST]
        BLASTDB=/opt/work/db/cdd_delta
        DATA_LOADERS=blastdb
        BLASTDB_PROT_DATA_LOADER=/opt/work/db/cdd_delta
        ;BLASTDB_NUCL_DATA_LOADER=
        
        ; Windowmasker settings
        ;[WINDOW_MASKER]
        ;WINDOW_MASKER_PATH=/opt/work/db/windowmasker
        ; end of file

    Run deltablast locally::

        deltablast -query 9mht.fasta -db cdd_delta

    **Example:**

    >>> bilab.webservice.NCBIStandalone('test',
                        queryid='test',
                        evalue=1e-10)
    output to file: default in xml
    
    >>> bilab.webservice.NCBIStandalone('test',
                        queryid='test',
                        filename='out',
                        gopen=10,
                        gext=2,
                        matrix='BLOSUM62',
                        outfmt=5,
                        pseudocount=0,
                        inclusion_ethresh=0.002,
                        domain_inclusion_ethresh=0.05,
                        num_iterations=5,
                        evalue=1e-10)
    
    Threads: using the number of logical cpus

    """
    if sequence == 'test':
        queryid='test'
        sequence=('MIEIKDKQLTGLRFIDLFAGLGGFRLALESCGAECVYSNEWDKYAQEVYEMNFGEKPEGD'
            'ITQVNEKTIPDHDILCAGFPCQAFSISGKQKGFEDSRGTLFFDIARIVREKKPKVVFMENVKNFAS'
            'HDNGNTLEVVKNTMNELDYSFHAKVLNALDYGIPQKRERIYMICFRNDLNIQNFQFPKPFELNTFV'
            'KDLLLPDSEVEHLVIDRKDLVMTNQEIEQTTPKTVRLGIVGKGGQGERIYSTRGIAITLSAYGGGI'
            'FAKTGGYLVNGKTRKLHPRECARVMGYPDSYKVHPSTSQAYKQFGNSVVINVLQYIAYNIGSSLNF'
            'KPY')
    try:
        sequence = ''.join(sequence.split())
        _ = sequence.isalpha()

    except AttributeError:
        raise ValueError('sequence must be a string.')
    else:
        if not _:
            raise ValueError('not a valid sequence')
        else:
            query_tmp_file = generate_temp_query(queryid, sequence)

    ncbi_cmd = []
    # basic parameters
    evalue = float(kwargs.pop('evalue', 1e-10))
    db = str(kwargs.pop('db','cdd_delta'))
    gapopen = int(kwargs.pop('gopen', 10))
    gapextend = int(kwargs.pop('gext', 2))
    matrix = str(kwargs.pop('matrix','BLOSUM62'))
    outfmt = int(kwargs.pop('outfmt', 5))
    if outfmt == 5:
        suffix = '.xml'
    elif outfmt == 3:
        suffix = '.txt'
    # Multiple threadings
    numThreads = get_num_logical_cpus()
    # PSSM engine options
    pseudocount = int(kwargs.pop('pseudocount', 0))
    pw_align_ethresh = float(kwargs.pop('inclusion_ethresh', 0.002))
    domain_inclusion_ethresh = float(kwargs.pop('domain_inclusion_ethresh', 0.05))
    # PSI-blast options
    num_iterations = int(kwargs.pop('num_iterations', 5))
    out_pssm = bool(kwargs.pop('out_pssm', True))
    out_ascii_pssm = bool(kwargs.pop('out_ascii_pssm', True))

    # Generate a command
    ncbi_cmd.append(program)
    ncbi_cmd.extend(['-query', str(query_tmp_file)])
    ncbi_cmd.extend(['-evalue',str(evalue)])
    ncbi_cmd.extend(['-db', str(db)])
    ncbi_cmd.extend(['-gapopen', str(gapopen)])
    ncbi_cmd.extend(['-gapextend', str(gapextend)])
    ncbi_cmd.extend(['-matrix', str(matrix)])
    ncbi_cmd.extend(['-outfmt', str(outfmt)])
    ncbi_cmd.extend(['-pseudocount', str(pseudocount)])
    ncbi_cmd.extend(['-inclusion_ethresh', str(pw_align_ethresh)])
    ncbi_cmd.extend(['-domain_inclusion_ethresh', str(domain_inclusion_ethresh)])
    ncbi_cmd.extend(['-num_iterations', str(num_iterations)])
    ncbi_cmd.extend(['-num_threads', str(numThreads)])

    if out_pssm:
        if filename is not None:
            out_pssm_file = filename + '_check_point.chk'
        else:
            out_pssm_file = program + '_check_point.chk'
        ncbi_cmd.extend(['-out_pssm', out_pssm_file])
    if out_ascii_pssm:
        if filename is not None:
            out_ascii_pssm_file = filename + '_pssm.txt'
        else:
            out_ascii_pssm_file = program + '_pssm.txt'
        ncbi_cmd.extend(['-out_ascii_pssm', out_ascii_pssm_file])
    if filename is not None:
        filename += suffix
        ncbi_cmd.extend(['-out', filename])
        # launch the command
        run(ncbi_cmd, verbose = True)
        # read xml file into pssm
        pssm_obj = PSSM(filename, sequence=sequence,
                        rawscorefile=out_ascii_pssm_file)
    else:
        # catch up from stdout
        output = run(ncbi_cmd, verbose = True)
        pssm_obj = PSSM(output, sequence=sequence,
                        rawscorefile=out_ascii_pssm_file)

    if deletetmp and os.path.exists(query_tmp_file):
        if verbose:
            print('Delete tmp query file:{0}'.format(query_tmp_file))
        os.unlink(query_tmp_file)

