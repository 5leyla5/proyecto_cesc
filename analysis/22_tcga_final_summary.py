import csv
from pathlib import Path

OUT = Path("results/tcga/tcga_final_project_summary.csv")

rows = [
    {
        "analysis": "Mutation burden",
        "comparison": "G2 vs G3",
        "group_1": "G2",
        "group_1_result": "median = 69 mutations",
        "group_2": "G3",
        "group_2_result": "median = 75 mutations",
        "p_value": "0.8079705",
        "fdr": "",
        "conclusion": "No significant difference"
    },
    {
        "analysis": "Mutation burden",
        "comparison": "Squamous vs Adenocarcinoma",
        "group_1": "Squamous",
        "group_1_result": "median = 74 mutations",
        "group_2": "Adenocarcinoma",
        "group_2_result": "median = 63 mutations",
        "p_value": "0.06856431",
        "fdr": "",
        "conclusion": "No significant difference; nominal trend"
    },
    {
        "analysis": "Recurrent genes",
        "comparison": "G2 vs G3",
        "group_1": "G2",
        "group_1_result": "Nominal gene-level differences observed",
        "group_2": "G3",
        "group_2_result": "Nominal gene-level differences observed",
        "p_value": "",
        "fdr": "All tested genes FDR >= 0.05",
        "conclusion": "No robust gene-level differences after multiple-testing correction"
    },
    {
        "analysis": "Recurrent genes",
        "comparison": "Squamous vs Adenocarcinoma",
        "group_1": "Squamous",
        "group_1_result": "KRAS 2.54%; TP53 4.24%; ERBB2 2.97%",
        "group_2": "Adenocarcinoma",
        "group_2_result": "KRAS 19.6%; TP53 17.4%; ERBB2 13.0%",
        "p_value": "Nominal signals",
        "fdr": "Not significant after FDR",
        "conclusion": "Exploratory differences only"
    },
    {
        "analysis": "SBS6 spectrum",
        "comparison": "Whole TCGA-CESC cohort",
        "group_1": "C>T",
        "group_1_result": "60.28%",
        "group_2": "C>G",
        "group_2_result": "20.59%",
        "p_value": "",
        "fdr": "",
        "conclusion": "C>T is the predominant substitution class"
    },
    {
        "analysis": "Ti/Tv",
        "comparison": "Whole TCGA-CESC cohort",
        "group_1": "Transitions",
        "group_1_result": "44,006",
        "group_2": "Transversions",
        "group_2_result": "22,640",
        "p_value": "",
        "fdr": "",
        "conclusion": "Ti/Tv = 1.9437"
    },
    {
        "analysis": "APOBEC signatures (SBS2 + SBS13)",
        "comparison": "G2 vs G3",
        "group_1": "G2",
        "group_1_result": "median = 47.02%",
        "group_2": "G3",
        "group_2_result": "median = 57.63%",
        "p_value": "0.03576317",
        "fdr": "0.3576317",
        "conclusion": "Higher descriptive APOBEC activity in G3, but not significant after FDR"
    },
    {
        "analysis": "APOBEC signatures (SBS2 + SBS13)",
        "comparison": "Adenocarcinoma vs Squamous",
        "group_1": "Adenocarcinoma",
        "group_1_result": "median = 47.50%",
        "group_2": "Squamous",
        "group_2_result": "median = 52.82%",
        "p_value": "0.6357358",
        "fdr": "0.8369486",
        "conclusion": "No significant histological difference"
    },
    {
        "analysis": "SBS40",
        "comparison": "Adenocarcinoma vs Squamous",
        "group_1": "Adenocarcinoma",
        "group_1_result": "median = 3.37%; SBS40 >=1% in 51.7%",
        "group_2": "Squamous",
        "group_2_result": "median = 0%; SBS40 >=1% in 21.6%",
        "p_value": "0.000763595",
        "fdr": "0.01527191",
        "conclusion": "Significantly greater SBS40 contribution in adenocarcinoma"
    },
    {
        "analysis": "SBS40 sensitivity",
        "comparison": "Adenocarcinoma vs Squamous; >=100 SNV and cosine >=0.90",
        "group_1": "Adenocarcinoma",
        "group_1_result": "9/20 SBS40 >0",
        "group_2": "Squamous",
        "group_2_result": "18/119 SBS40 >0",
        "p_value": "Wilcoxon = 0.000820632",
        "fdr": "",
        "conclusion": "SBS40 association remains in stricter sensitivity analysis"
    }
]

OUT.parent.mkdir(parents=True, exist_ok=True)

with OUT.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "analysis",
            "comparison",
            "group_1",
            "group_1_result",
            "group_2",
            "group_2_result",
            "p_value",
            "fdr",
            "conclusion"
        ]
    )
    writer.writeheader()
    writer.writerows(rows)

print("\n===== RESUMEN FINAL TCGA-CESC =====\n")

for r in rows:
    print(f"{r['analysis']}")
    print(f"  {r['comparison']}")
    print(f"  -> {r['conclusion']}\n")

print("Archivo creado:")
print(OUT)
