def getLPS(pattern):
    """
    Computes the longest proper prefix which is also a suffix (LPS) for each prefix of the pattern.

    Args:
        pattern (str): The pattern string for which LPS needs to be computed.

    Returns:
        A list of integers representing the LPS for each prefix of the pattern.
    """
    len_pat = len(pattern)
    lps = [0] * len_pat
    j = 0
    for i in range(1, len_pat):
        while j > 0 and pattern[i] != pattern[j]:
            j = lps[j-1]
        if pattern[i] == pattern[j]:
            j += 1
            lps[i] = j
    return lps

def kmp(text, pattern):
    """
    Finds all occurrences of the pattern in the text using the Knuth-Morris-Pratt algorithm.

    Args:
        text (str): The input text in which pattern needs to be searched.
        pattern (str): The pattern string that needs to be searched in the input string.

    Returns:
        Returns pattern and indices of the pattern if it found. None otherwise.
    """

    
    len_text, len_pat = len(text), len(pattern)
    lps = getLPS(pattern)
    i, j = 0, 0
    locs = []
    while i < len_text:
        if pattern[j] == text[i]:
            i += 1
            j += 1
            if j == len_pat:
                locs.append(i - j)
                j = lps[j-1]
        else:
            if j > 0:
                j = lps[j-1]
            else:
                i += 1
    return (pattern, locs) if locs else None
   

def write_result(pattern, locs, file, lock):
    """
    Writes the pattern to file.

    Args:
        pattern (str): The found pattern string.
        lock (multiprocessing.Lock): A lock object for synchronizing writes.

    Returns:
        None.
    """

    lock.acquire()
    file.write(pattern + '\t' + str(len(locs)) + '\t' + str(locs) + '\n')
    lock.release()

def search(pattern_file, text_file, output_file, ignore_case):
    """
    Searches for all patterns in the given text file.

    Args:
        pattern_file (str): The name of the file containing the patterns to be searched.
        text_file (str): The name of the file containing the text in which patterns needs to be searched.
        output_file (str): The name/directory of the output file.
        ignore_case (bool): A flag indicating whether to ignore case sensitivity while searching for the patterns.

    Returns:
        None.
    """
    
    #Check that both input files are fasta files.
    if not pattern_file.endswith((".fasta")):
        raise ValueError("pattern_file must be in .fasta format")   
    if not text_file.endswith((".fasta")):
        raise ValueError("text_file must be in .fasta format")   
    
    #Read fasta files with SeqIO and parse into patterns, removing capitalization if inputed.
    with open(pattern_file) as f:
        if ignore_case:
            patterns = [str(record.seq.lower()) for record in SeqIO.parse(f, "fasta")]
        else:
            patterns = [str(record.seq) for record in SeqIO.parse(f, "fasta")]

    with open(text_file) as f:
        genome = SeqIO.read(f,"fasta")
        print(genome)
        if ignore_case:
            genome = genome.seq.lower()
        else:
            genome = genome.seq

    #Open the output file and run the code, each pattern with its own process.
    with open(output_file, "w+") if output_file else sys.stdout as f:
        lock = multiprocessing.Lock()
        pool = multiprocessing.Pool()
        f.write('Sequences found in genome: \n' +
                'SEQ \t COUNT \t LOCS \n')
        for pattern in patterns:
            pool.apply_async(kmp, (genome, pattern), callback=lambda x: (write_result(x[0], x[1], f, lock) if x else None))
        pool.close()
        pool.join()


if __name__ == '__main__':
    import argparse, sys, multiprocessing
    from Bio import SeqIO

    parser = argparse.ArgumentParser(description='The KMP tool searches for a short sequences in a long text. '\
                                      'The tool outputs all the sequences that were found in the text.')
    parser.add_argument('pattern_file', help = 'Fasta file containing a set of short sequences')
    parser.add_argument('text_file', help = 'Fasta file containing text to search in.')
    parser.add_argument('-o', '--output_file', default = None, help = 'Optional file to write search results to (default: outputs to command line).')
    parser.add_argument('--ignore-case', action = 'store_true', help = 'Ignore case when searching')
    args = parser.parse_args()

    # if not os.path.isfile(args.pattern_file):
    #     print(f"Error: {args.pattern_file}: No such file or directory.", file=sys.stderr)
    #     sys.exit(1)

    # if not os.path.isfile(args.text_file):
    #     print(f"Error: {args.text_file}: No such file or directory.", file=sys.stderr)
    #     sys.exit(1)

    search(args.pattern_file, args.text_file, args.output_file, args.ignore_case)
