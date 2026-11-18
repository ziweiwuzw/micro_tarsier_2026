# Workflow order

The numbered folders follow the Materials and Methods. Each entry point submits or runs the scripts beside it.

| Order | Step | Entry point | Output |
| --- | --- | --- | --- |
| 1.1 | Clean paired reads | `01_reads_qc_host/process_reads_fastp.sh` | Clean FASTQ files and fastp reports |
| 1.2 | Remove host reads | `01_reads_qc_host/remove_host_reads.slurm` | Paired nonhost FASTQ files |
| 1.3 | Summarize samples | `01_reads_qc_host/summarize_fastp.sh`, `summarize_nonhost.slurm`, `build_dataset_summary.py` | Sample summary tables |
| 1.4 | Classify reads from four focal samples | `01_reads_qc_host/ccmetagen/` | CCMetagen assignments from fastp-cleaned reads |
| 2.1 | Assemble each sample | `02_assembly_binning/assemble_samples.slurm` | MEGAHIT contigs |
| 2.2 | Match assemblies to nonhost reads | `02_assembly_binning/prepare_sample_manifest.py` | Shared sample manifest |
| 2.3a | Bin host-group contigs with SemiBin2 and VAMB | `02_assembly_binning/submit_joint_binning.sh` | Candidate bins |
| 2.3b | Bin individual assemblies with COMEBin | `02_assembly_binning/submit_comebin.sh` | Candidate bins from same-species Strobealign coverage |
| 2.4 | Refine bins | `02_assembly_binning/refinement/submit.sh` | Final candidate MAGs |
| 3 | Assess MAGs, form SGBs and estimate abundance | `03_mag_sgb_kma/submit.sh check`, then `submit.sh submit` | MAG quality, SGB IDs, taxonomy and KMA abundance |
| 4 | Annotate SGB member MAGs | `04_sgb_annotation/submit.sh all`; `submit.sh fegenie` for iron genes | Gene calls and functional tables |
| 5 | Analyse host gene families | `05_host_cafe5/prepare_cafe5_inputs.py`, `run_cafe5.slurm`, `extract_focal_changes.py` | CAFE5 family and branch results |

The two binning paths in step 2.3 run from the same prepared samples and meet at bin refinement. KMA profiling and functional annotation both start after SGB IDs are assigned. The CCMetagen branch uses fastp-cleaned reads without host read removal.
