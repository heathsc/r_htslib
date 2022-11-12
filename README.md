# r_htslib
Lightweight wrapper around htslib for access to SAM/BAM/CRAM/BCF/VCF files.

Changes

 - 0.12.0. Remove the mmalloc function as they are superseded and safer options are available
 - 0.11.0. Add many new functions in vcf from gemBS-rs.  Also breaking change to Hts::open where now the file parameter is an Option<>
 - 0.10.2. Add several accessor functions to Hts:: and HtsHdr:: so that the underlying SAM/VCF headers do not need to be directly accessed
 - 0.10.1. Implement Send + Sync for faidx::Sequence
 - 0.10.0. Change faidx so that (a) sub chromosome regions can be fetched and (b) uses PhantomData to flag that Sequence owns data  
 - 0.9.3. Add ability to specify . or * for regions in Hts::make_region_list()
 - 0.9.2. Add seq_lengths() function to hts that obtains a vector of sequence lengths from sam/bam/cram/vcf/bcf files, and returns an empty vector otherwise
