import polars as pl
import pathlib

"""These are all the iron uptake and transport genes found in literature associated with Cefiderocol resistance """

genes = ["fhuA", "fepA", "fbpA", "efeO", "exbB", "exbD", "fiuA", "fur", "iutA",
         "baeS", "envZ", "cirA", "feoA", "sitC", "apbC", "fepG", "fepC", "fetB",
         "fetA", "fecA", "fiu", "ompR", "yicI", "yicJ", "yicL", "chrA", "OmpR", "tonB", "sugE", "yicM", "PmrB", "ompK35", "ompK36", "ompK37"]

snippy_folder = pathlib.Path("/home/raulc/Desktop/Test_Cefiderocol/out/snippy_kp48it")
cv2_filepath = snippy_folder / "CV2/snps.tab"
cv2 = pl.read_csv(cv2_filepath, separator="\t")

iron_genes = []

for i in genes:
    df2 = cv2.filter(pl.col("GENE") == i)
    if df2.shape[0] != 0:
        iron_genes.append(df2)
if len(iron_genes) != 0:
    scv2 = pl.concat(iron_genes)
    fn = "iron_genes_cv2.csv"
    filepath = snippy_folder / fn
    scv2.write_csv(filepath, separator="\t")
else:
    print(f"No mutations found in iron uptake and transport genes for CV2!")


cv4_filepath = snippy_folder / "CV4/snps.tab"
cv4 = pl.read_csv(cv4_filepath, separator="\t")

iron_genes = []

for i in genes:
    df4 = cv4.filter(pl.col("GENE") == i)
    if df4.shape[0] != 0:
        iron_genes.append(df4)
if len(iron_genes) != 0:
    scv4 = pl.concat(iron_genes)
    fn = "iron_genes_cv4.csv"
    filepath = snippy_folder / fn
    scv4.write_csv(filepath, separator="\t")
else:
    print(f"No mutations found in iron uptake and transport genes for CV4!")

cv5_filepath = snippy_folder / "CV5/snps.tab"
cv5 = pl.read_csv(cv5_filepath, separator="\t")

iron_genes = []

for i in genes:
    df5 = cv5.filter(pl.col("GENE") == i)
    if df5.shape[0] != 0:
        iron_genes.append(df5)
if len(iron_genes) != 0:
    scv5 = pl.concat(iron_genes)
    fn = "iron_genes_cv5.csv"
    filepath = path / fn
    scv5.write_csv(filepath, separator="\t")
else:
    print(f"No mutations found in iron uptake and transport genes for CV5!")