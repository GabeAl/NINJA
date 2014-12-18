NINJA
=====

NINJA Is Not Just Another Reference Sequence Matcher

This is a new ultra-high-speed OTU picking pipeline. 

Here's the rundown. 

1. Make reference DB by providing a set of OTUs (i.e. f/greengenes) and optional taxa file (also f/greengenes). This is done with ninja_prep.

2. Run the resulting output fasta through bowtie2's indexer, generating the bw index.
    bowtie2-align-s --no-head --no-unal -o3 -p4 -f reads.fna -x Ninja97 -S align.txt --mp "1,1"--rdg "1,1" --rfg "1,1" --score-min "L,0,-.03" -k 1 --norc --fast

3. Run ninja_filter on your reads (or don't -- but filtering makes things faster and cleaner).

4. Run your reads through bowtie against the index generated in step 2.

5. If you did step 3, run ninja_parse on bowtie2's sam output. Otherwise run ninja_parse_s on it.

6. ???

7. profit 

