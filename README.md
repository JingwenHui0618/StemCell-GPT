# StemCell-GPT Pipeline

This repository contains a skeleton for a CRISPR guide design pipeline targeting human stem cells. The code here demonstrates how independent modules can be orchestrated into a single workflow.

The existing modules referenced in the code are not included in this repository. Replace the stubs with the actual implementations from your private codebase.

## Workflow Overview
1. **User Input** – Choose cell type (HSPC or WT), specify target region/gene, and optionally supply a personal VCF file.
2. **Guide Collection** – Use `guide_filter.py` to generate candidate sgRNAs while avoiding common SNPs. If a personal VCF is provided, run `personal_snp_filter.py` to remove guides overlapping personal variants.
3. **Off-target Scoring** – Score each guide using an off-target tool (e.g., the CRISPOR API). Select up to 20 guides with the best off-target profiles.
4. **On-target Scoring** – Compute a raw CRISPR score for each guide, embed the sequence using DNABERT, and evaluate the fine‑tuned on-target model stored in `OT_model/`.
5. **Forecast Analysis** – Predict the indel distribution with Forecast and calculate an MMEJ/NHEJ score.
6. **Final Ranking** – Combine the fine‑tuned on‑target score and Forecast score using `0.3 * on_target + 0.7 * forecast` and rank guides.
7. **ssODN Design** – Use the `ssodn_design.py` module to design a repair template for the top guide(s).
8. **Result Output** – Present ranked guides and the designed ssODN to the user.

## Usage
Set the environment variable `OPENAI_API_KEY` with your OpenAI key before running the pipeline:
```bash
export OPENAI_API_KEY=<your-key-here>
```
Then run:
```bash
python pipeline.py --cell-type HSPC --region chr1:1000-2000 --gene TP53 --vcf my.vcf
```

## Disclaimer
The repository lacks the large genomic data files and the actual implementation of the individual modules. The provided `pipeline.py` serves as an example of how those components can be tied together.
