# Gene Finding NHEMMO (Non Homogeneous Empirical Markov Model for Organisms)

## ABSTRACT
Our project is a re-implementation of GlimmerHMM, a gene finding software. It is designed to predict gene enconding regions in genetic data, specifically whether a portion of an organism's genome represents an exon (initial, internal, final, single), intron (phase 0, phase 1, phase 2), or intergenic region.

Our model is trained on genomic data from Clown anemonefish (_Amphiprion ocellaris_) chromosomes, in light of our _Finding Nemo_ theme, and can be used to either infer the region-based genomic breakdown of other chromosomes of this species, or that of similiar species via transfer learning. Users may also train their own model using the training pipeline detailed below.

## Installation
1. Clone this repository: `git clone https://github.coecis.cornell.edu/ahc248/gene_finding_nhemmo.git`
2. Install the following dependencies:
   1. `pip install numpy`
   2. `pip install python-dotenv`

## Usage

### Inference via Provided Clown anemonefish (_Amphiprion ocellaris_) Model:
1. Run `python region_inference.py`.
2. Locate results in the `results` directory.

**NOTE: All needed training files are already located in their correct locations, listed below:**
   - Exon file: `fasta_gtf_exon_files` directory
   - Transition & emission matrices: `matrices` directory


### Downloading Clown anemonefish (_Amphiprion ocellaris_) Data
**NOTE: `.gtf` and `.fa` files are too large to store in this repository, so download instructions are provided below.**
1. Download files from [_Ensembl_](https://useast.ensembl.org/Amphiprion_ocellaris/Info/Index)
   1. **Download the FASTA file:**
      1. Under _Gene annotation_ click [`Download FASTA`](https://ftp.ensembl.org/pub/release-113/fasta/amphiprion_ocellaris/).
      2. Click [`dna/`](https://ftp.ensembl.org/pub/release-113/fasta/amphiprion_ocellaris/dna/) folder (NOTE: the heading of the page should now be `Index of /pub/release-113/fasta/amphiprion_ocellaris`).
      3. Click `Amphiprion_ocellaris.ASM2253959v1.dna.toplevel.fa` to download the .zip file.
      4. Extract the .zip file into the `fasta_gtf_exon_files` directory.
      5. Remove `Amphiprion_ocellaris.ASM2253959v1.dna.toplevel.fa` from its parent directory (it should not be located within a subdirectory, but rather directly inside `fasta_gtf_exon_files`). 
   2. **Download the GTF file:**
      1. Under _Gene annotation_ click [`Download GTF`](https://ftp.ensembl.org/pub/release-113/gtf/amphiprion_ocellaris/).
      2. Click `Amphiprion_ocellaris.ASM2253959v1.113.gtf` to download the .zip file (NOTE: the heading of the page should now be `Index of /pub/release-113/gtf/amphiprion_ocellaris`).
      3. Extract the .zip file into the `fasta_gtf_exon_files` directory.
      4. Remove `Amphiprion_ocellaris.ASM2253959v1.113.gtf` from its parent directory (it should not be located within a subdirectory, but rather directly inside `fasta_gtf_exon_files`). 

### Full Pipeline: Training & Inference via Provided Data
1. Add FASTA (`.fa`) and GTF (`.gtf`) files for organism of interest to `fasta_gtf_exon_files` directory.
    - Example FASTA filename: `Amphiprion_ocellaris.ASM2253959v1.dna.toplevel.fa`
    - Example GTF filename: `Amphiprion_ocellaris.ASM2253959v1.113.gtf`
2. Update environment (`.env`) variables `FASTA_FILEPATH`, `GTF_FILEPATH`, `OUTPUT_EXON_FILEPATH`, and `SUMMARY_FILEPATH` to specify the filepaths of the FASTA file, GTF file, resulting exon file, and resulting summary file, respectively.
3. Run `python gtf_to_exons.py` to generate the exon file and corresponding summary file.
4. Update environment (`.env`) variable `NUM_TRAINING_LABELS` to specify the desired Markov orders upon which to train the new model.
5. Delete all pre-existing `.npy` files in the `matrices` directory.
6. Run `python training.py` to train the new model and generate the corresponding transition and emission matrices (stored as `.npy` files).
7. Update environment (`.env`) variables `MAX_MARKOV_ORDER` and `INFERENCE_FILEPATH` to specify the number of chromosomes to allocate for the training set (the remaining chromosomes can be used for inference), and the filepath of the resulting results file, respectively.
8. Delete all pre-existing `.csv` files in the `results` directory.
9. Run `python region_inference` to infer the gene encoding regions in the organism's genome, using the new model.
10. Locate results in the `results` directory.
11. Run `analysis.py` to reproduce bar charts analyzing the accuracy of the infered gene regions.

Authors: Chris Gunning (clg238), April Collamore (ahc248), Maya Deen (mad362), Natalie Suggs (nrs86)
