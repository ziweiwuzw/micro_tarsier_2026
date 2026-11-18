# Primate metagenome analysis scripts

These scripts accompany the Materials and Methods of our primate metagenome study. The folders follow the order of the analysis. [Table 1](Table_1_software_parameters.md) lists the software versions and main settings.

| Step | Analysis | Main entry point |
| --- | --- | --- |
| 1 | Read cleaning, host read removal and sample summaries | `01_reads_qc_host/process_reads_fastp.sh`, then `remove_host_reads.slurm` |
| 2 | Assembly, coverage mapping, binning and bin refinement | `02_assembly_binning/assemble_samples.slurm`, `submit_joint_binning.sh`, `submit_comebin.sh`, then `refinement/submit.sh` |
| 3 | MAG quality, SGB construction, taxonomy and abundance | `03_mag_sgb_kma/submit.sh` |
| 4 | Gene prediction and functional annotation | `04_sgb_annotation/submit.sh` |
| 5 | Host gene family analysis | `05_host_cafe5/prepare_cafe5_inputs.py`, then `run_cafe5.slurm` |

[Workflow order](WORKFLOW_ORDER.md) shows the branches and the inputs passed between steps. Related commands share an entry point where they are run together; supporting scripts remain beside that entry point.

## Running the scripts

These scripts were written for a Slurm cluster. Paths marked `/path/to/` are placeholders. Set them to your input, reference database, Conda and output locations before running the scripts. Set `SBATCH_ACCOUNT` and the Slurm partition for your cluster. Raw reads, reference databases and analysis outputs are not included.

For example, set `PROJECT_DIR` and `CONDA_SH` before `01_reads_qc_host/process_reads_fastp.sh`. To build the assembly manifest, pass `--assembly-root`, one or more `--read-root` values, and `--output` to `02_assembly_binning/prepare_sample_manifest.py`.

The four focal samples also have a CCMetagen branch in `01_reads_qc_host/ccmetagen/`. That branch uses fastp-cleaned reads without host read removal. For COMEBin, `submit_comebin.sh` runs CPU jobs with 64 threads and eight workers; missing same-species BAMs are made with Strobealign.
