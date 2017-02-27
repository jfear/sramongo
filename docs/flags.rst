.. _database_flags:

============================================================
Database Flags
============================================================

``sra2mongo`` adds several flags based on data found in the SRA metadata.  These
flags are created to make querying specific metadata easier. Flags can be at
either the `Ncbi`_ level or the `Run`_ level.

Ncbi
====

SE
    Single-end reads according to the SRA ``library_layout``.

PE
    Paired-end reads according to the SRA ``library_layout``.

RNASeq
    SRA gives evidence that the experiment should be RNA-seq. The flag is added
    if ``library_source`` is ``TRANSCRIPTOMIC`` or if ``library_strategy`` is
    ``RNA-Seq``.


Run
===

SE
    Single-end reads according to the SRA.

PE
    Paired-end reads according to the SRA.

no_read_information
    nreads in SRA was "variable" and no read information was provided.

pe_reads_not_equal_len
    Read 1 and read 2 lengths are more than 5bp different according to the SRA.

pe_reads_not_equal_count
    Read 1 and read 2 counts are more than 10 reads different according to the
    SRA.
