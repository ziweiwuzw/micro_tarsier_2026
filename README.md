# Primate metagenome Methods: curated server scripts

This draft groups the scripts by the order of the manuscript Methods. Only stage directories are numbered; script names describe their function, and [WORKFLOW_ORDER.md](WORKFLOW_ORDER.md) gives the execution and branch order. `SOURCE_MANIFEST.tsv` records the original source and SHA-256 for each retained source file; `SOURCE_AUDIT.tsv` records excluded and merged copies. No original source directory was modified. Cluster data, databases, MAGs and execution logs are not included.

## Stage order and entry points

| Order | Stage | Entry point and outcome |
| --- | --- | --- |
| 1 | `01_reads_qc_host/` | `process_reads_fastp.sh` → `remove_host_reads.slurm`; summary scripts produce sample metadata and paired nonhost FASTQ. `ccmetagen/` is a separate four-sample branch using fastp-cleaned reads **without host removal**. |
| 2a | `02_assembly_binning/` | `assemble_samples.slurm` produces per-sample contigs ≥1 kb. |
| 2b | `02_assembly_binning/` | `prepare_sample_manifest.py` discovers assemblies and nonhost reads → `prepare_host_group_contigs.slurm` → `prepare_joint_binning.py` → `submit_joint_binning.sh` for host-group SemiBin2/VAMB bins. |
| 2c | `02_assembly_binning/` | `submit_comebin.sh` coordinates CPU mapping and GPU COMEBin binning for each sample. |
| 2d | `02_assembly_binning/refinement/` | `submit.sh` prepares groups and submits `run_metawrap_refinement.slurm`, producing the final candidate-bin source for stage 3. |
| 3 | `03_mag_sgb_kma/` | `submit.sh check` then `submit.sh submit`: CheckM2 → HQ MAGs → dRep/SGB IDs; taxonomy, source metadata, KMA abundance and annotation branch after SGB IDs. `stage_job.slurm` and `kma_job.slurm` each run several short modes. |
| 4 | `04_sgb_annotation/` | Stage 3 launches `launch_annotation.slurm` → `submit.sh all` for Prodigal, DRAM2, dbCAN, VFDB/Toxinome and source tables. `submit.sh fegenie` launches the separate FeGenie branch. `summary_job.slurm` combines short summary/export jobs. |
| 5 | `05_host_cafe5/` | `prepare_cafe5_inputs.py` → `run_cafe5.slurm` → `extract_focal_changes.py`; `prepare_tree_statistics.sh` and `plot_cafe5_tree.R` are distinct tree/stat and plotting helpers. |
| Optional | `optional/hecatomb_pilot/` | Two-host Hecatomb pilot; its contribution to the manuscript remains unverified. |

KMA and annotation depend on SGB IDs but can then run independently. PPanGGOLiN was excluded because the author confirmed it was not used.

## Running on Spartan

1. Review all absolute `/data/...`, `/home/...` and `/apps/...` paths, Slurm account/partitions, environment variables, databases and input manifests. Scripts encode the original cluster layout.
2. Run stage 1 and the stage 2 chains. `prepare_sample_manifest.py --assembly-root <assemblies> --read-root <reads> --output <manifest>` replaces the old global-catalogue preparation for the selected path. `--assembly-list` remains available for an existing list. COMEBin requires an explicit `MAPPER` value; the executed mapper still needs confirmation from run records.
3. Run `02_assembly_binning/refinement/submit.sh` with a fresh `PIPELINE_ROOT`. Then set `SOURCE_ROOT`, `INPUT_READS_MANIFEST`, `PIPELINE_ROOT`, `ANN_OUT` and site-specific values for stage 3. Check inputs with `bash 03_mag_sgb_kma/submit.sh check`, then submit with `bash 03_mag_sgb_kma/submit.sh submit`.
4. After ORF calls, run `bash 04_sgb_annotation/submit.sh fegenie` if FeGenie belongs to the final analysis.

Stage 3 requires 32 complete group directories and a new output root. The stage 2 refiner can pass through a single binner's bins, including a default *Papio ursinus* exception. Read each group's `.metawrap_refinement.done` before describing the collection as universally MetaWRAP-refined.

The publication-facing [Supplementary Table S1](docs/Supplementary_Table_S1_software.md) records software and settings; `PARAMETERS.md` points to that single table. `MERGE_AUDIT.md` explains which steps were consolidated, and `RECONCILIATION.md` lists claims requiring run evidence. Static checks pass locally; biological outputs and Slurm resource behavior need validation on Spartan.
