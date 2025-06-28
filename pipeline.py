import argparse
import json
import os
from typing import List, Dict

# Placeholder imports for the modules referenced in the design. Replace these with
# actual implementations from your private repository.
try:
    import guide_filter
    import personal_snp_filter
    import off_target
    import on_target
    import forecast
    import ssodn_design
except ImportError:
    # Stubs so the file remains runnable even without the real modules.
    class guide_filter:
        @staticmethod
        def design_guides(cell_type: str, region: str, gene: str) -> List[str]:
            """Return a list of candidate guide sequences."""
            return []

    class personal_snp_filter:
        @staticmethod
        def filter_guides(guides: List[str], vcf_file: str) -> List[str]:
            return guides

    class off_target:
        @staticmethod
        def score(guide: str) -> float:
            return 0.0

    class on_target:
        @staticmethod
        def raw_score(guide: str) -> float:
            return 0.0

        @staticmethod
        def finetune_score(guide: str, raw_score: float) -> float:
            return raw_score

    class forecast:
        @staticmethod
        def score(guide: str) -> float:
            return 0.0

    class ssodn_design:
        @staticmethod
        def design(guide: str, mutation: str) -> str:
            return ""


FINAL_SCORE_RATIO = 0.3
FORECAST_RATIO = 0.7

def rank_guides(guides: List[str]) -> List[Dict]:
    results = []
    for g in guides:
        ot_score = off_target.score(g)
        raw = on_target.raw_score(g)
        ft = on_target.finetune_score(g, raw)
        fc = forecast.score(g)
        final = FINAL_SCORE_RATIO * ft + FORECAST_RATIO * fc
        results.append({
            "guide": g,
            "off_target": ot_score,
            "raw_on_target": raw,
            "fine_tuned_on_target": ft,
            "forecast": fc,
            "final_score": final,
        })
    results.sort(key=lambda x: x["final_score"], reverse=True)
    return results

def main():
    parser = argparse.ArgumentParser(description="StemCell-GPT sgRNA pipeline")
    parser.add_argument("--cell-type", choices=["HSPC", "WT"], required=True)
    parser.add_argument("--region", required=True, help="Target genomic region")
    parser.add_argument("--gene", required=True, help="Gene name")
    parser.add_argument("--vcf", help="Optional personal VCF file")
    parser.add_argument("--mutation", default="", help="Description of desired edit")
    args = parser.parse_args()

    guides = guide_filter.design_guides(args.cell_type, args.region, args.gene)
    if args.vcf:
        guides = personal_snp_filter.filter_guides(guides, args.vcf)

    scored = rank_guides(guides)[:20]
    if scored:
        best = scored[0]["guide"]
        ssodn = ssodn_design.design(best, args.mutation)
    else:
        ssodn = ""

    output = {
        "guides": scored,
        "ssodn": ssodn,
    }
    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
