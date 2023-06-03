# KMP

In string computation, the exact pattern matching problem is the problem of
finding all the occurences of a pattern (string) P, in a text (string) S, where usually P is much
shorter than S. For example the pattern could be the world “stella” and the text the whole Divina
Commedia, or P can be the CCATTGTG motif and the text the human genome.
One strategy to speed up the computation is to create an index on the pattern P and use this index to
scan the text S in a more efficient way.
The Knuth-Morris-Pratt algorithm uses this approach. It first of all builds an index on P and then
uses it to scan S, applying simple rules to the index to decide how to shift the pattern.
This is a simple approach to apply this algorithm for scanning a genome fasta file and finding short sequences that appear in the genome. The output of this function is a file containing a list of sequences, their locations in the genome and their count.

# Installation
Make sure you have all required libraries
```bash
$ pip install -r requirements.txt
```

# Usage
The function has 2 mandatory inputs and 2 optional inputs:
- pattern_file: Fasta file containing a set of short sequences.
- text_file: Fasta file containing text to search in.
- -o, --output_file: Optional file to write search results to (default: outputs to command line)
- --ignore-case: Ignore case when searching

```bash
$ python3 KMP.py pattern_file text_file -o out.txt --ignore_case
```



