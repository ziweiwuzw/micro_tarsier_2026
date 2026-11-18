# Table 1. Software and main parameters

This table lets readers find the main settings without opening the scripts. Thread counts are per job.

| Stage | Software (version) | Purpose | Main settings |
| --- | --- | --- | --- |
| Read cleaning | fastp (0.22.0) | Remove adapters and low-quality reads | Paired-end adapter detection; Phred ≥25; ≤40% low-quality bases; ≤3 ambiguous bases; read length ≥36 bp; 12 threads. |
| Host removal | BWA-MEM (0.7.17) | Align reads to the host reference | `bwa mem -t 24`. |
| Host removal | SAMtools (1.21) | Keep pairs with neither mate mapped | `collate -u -O`; `fastq -f 12 -F 0x900 -n`; 8 threads. |
| Focal sample taxonomy | KMA  | Classify fastp-cleaned reads against NCBI nt | KMA `-ipe -1t1 -mem_mode -and -apm p -ef`; 48 threads. |
| Assembly | MEGAHIT (1.2.9) | Assemble each sample | Paired reads (`-1`, `-2`); `--min-contig-len 1000`; 32 threads. |
| Coverage | Strobealign (0.17.0) | Map reads to assembled contigs | Read-length profile (`-r`); `-R 6`; profile-specific index. Host-group jobs used 24 mapping threads; COMEBin mapping used 56. |
| Coverage | SAMtools (1.24) | Produce sorted BAM files | Exclude unmapped reads (`-F 4`); coordinate sort and index. Host-group jobs used 4 threads; COMEBin sorting used 8. |
| Host-group binning | SemiBin2 (2.4.1) | Recover bins from host-group contigs | `multi_easy_bin --self-supervised`. |
| Host-group binning | VAMB (5.0.4) | Recover bins from host-group contigs | `vamb bin default`. |
| Per-sample binning | COMEBin (1.1.0) | Bin each sample using same-species coverage | Contigs ≥2 kb; Strobealign BAMs; `-n 6 -d cpu -t 64 -s 42 -w 8`. |
| Bin refinement | MetaWRAP (1.3) | Refine candidate bins | `bin_refinement -c 50 -x 10 --quick`; 24 threads; 180 GB tool memory. |
| Bin refinement | CheckM (1.0.18) | Assess bins within MetaWRAP | MetaWRAP-compatible CheckM database. |
| MAG quality | CheckM2 (1.0.1) | Estimate completeness and contamination | Input extension `.fa`; 24 threads. High-quality MAGs: completeness ≥90%, contamination ≤5%. |
| SGB construction | dRep (3.6.2) | Dereplicate MAGs | Primary ANI 90%; secondary ANI 95%; minimum aligned fraction 0.30; fastANI; 24 threads. |
| SGB taxonomy | GTDB-Tk (2.4.0) | Classify SGB representatives | GTDB R220; `classify_wf`; Mash database; 32 CPUs and 16 pplacer CPUs. |
| SGB abundance | KMA (1.6.18) | Map nonhost read pairs to SGB representatives | `-ipe -t_db -mem_mode -ef -1t1 -nf -nc`; 16 threads per sample. SGB TPM uses summed contig `fragmentCount` and representative genome length. |
| Gene prediction | Prodigal (2.6.3) | Predict coding genes in member MAGs | `-p meta -g 11`; protein, nucleotide and GFF outputs. |
| Functional annotation | DRAM2 (2.0.0-beta38) | Annotate genes and metabolic functions | All member MAGs; `--input_genes --annotate --summarize --sum_ecos gut`; `dram_db,kofam,merops,metals,sulfur,card`. KOfam uses KO-specific score thresholds. |
| DRAM2 execution | Nextflow (25.10.4) | Run the DRAM2 workflow | Apptainer and Slurm profile; DRAM2 revision `4536cd93c9`. |
| DRAM2 execution | Apptainer (1.3.3) | Run DRAM2 containers | Used through the Nextflow Apptainer profile. |
| CAZyme annotation | dbCAN (5.2.9) | Assign CAZyme families | Protein mode; DIAMOND, HMM and dbCANsub; 8 threads per MAG. A family enters the selected SGB table when two methods support it. |
| Virulence and toxin screening | DIAMOND (2.2.8) | Search VFDB and Toxinome | VFDB: `blastp --very-sensitive`, identity and both coverages ≥70%, E ≤10⁻¹⁰. Toxinome: `blastp --sensitive`, identity ≥40%, query coverage ≥90%, subject coverage ≥60%, E ≤0.001. Both: ≤25 targets and 8 threads per MAG. |
| Iron metabolism | FeGenie (version not recorded) | Find iron uptake, regulation, oxidation, reduction and storage genes | Metagenome mode (`--meta`); 4 threads per MAG job. |
| Host gene families | CAFE5 (version not recorded) | Test expansion and contraction of primate gene families | Gamma model with three rate categories; significant branch changes at P ≤0.05. |

FeGenie and CAFE5 version numbers were not recorded in the supplied parameter list.
