# %%
import time
section_start = time.time()

from typing_extensions import ParamSpecArgs
import os
import sys
import argparse
import cProfile
import numpy as np
import pandas as pd
import scanpy as sc

print(f"[PROFILE] imports: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

parser = argparse.ArgumentParser(description="Process PBMC datasets with profiling.")
parser.add_argument("--data-dir", type=str, default="data",
                    help="Directory containing dataset subdirectories")
parser.add_argument("--data-set", type=str, default="pbmc3k",
                    help="Dataset name/subdirectory")
parser.add_argument("--out-dir", type=str, default="data",
                    help="Output directory")
parser.add_argument("--num-threads", type=int, default=1,
                    help="Number of threads")

args = parser.parse_args()

datadir = args.data_dir if args.data_dir.endswith("/") else args.data_dir + "/"
dataset = args.data_set
outdir = args.out_dir if args.out_dir.endswith("/") else args.out_dir + "/"
nthreads = args.num_threads

os.makedirs(outdir, exist_ok=True)

print(f"[PROFILE] argument_parsing: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")
sc.settings.n_jobs = nthreads

print(f"Dataset: {dataset}")
print(f"Using {sc.settings.n_jobs} Scanpy jobs")
print(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', 'not set')}")
print(f"MKL_NUM_THREADS={os.environ.get('MKL_NUM_THREADS', 'not set')}")
print(f"[PROFILE] scanpy_settings: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

# I/O
results_file = os.path.join(outdir, dataset + ".scanpy.h5ad")
profile_file = os.path.join(outdir, dataset + ".rank_genes_groups.prof")

adata = sc.read_10x_mtx(
    os.path.join(datadir, dataset, "filtered_gene_bc_matrices"),
    var_names="gene_symbols",
    cache=True,
)
adata.var_names_make_unique()

print(f"[PROFILE] io_read_10x: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

# preprocessing: basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"[PROFILE] basic_filtering: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()
#%%
# metric
#adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
#sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# filtering by slicing the AnnData object
#adata = adata[adata.obs.n_genes_by_counts < 2500, :]
#adata = adata[adata.obs.pct_counts_mt < 5, :]


# normalize to 10K reads per cell
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print(f"[PROFILE] normalization_log1p: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

# highly variable genes
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# freeze data
adata.raw = adata

# retain highly variable genes
adata = adata[:, adata.var.highly_variable]

print(f"[PROFILE] highly_variable_genes: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()
# regres out effects of total counts per cell an d% mitochondrial genes
#sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata)

print(f"[PROFILE] scaling: {time.time() - section_start:.6f} s")

# %%
# report adata - so we can check ot see if we are comparable to Seurat
# adata.write(results_file)
# adata

# %%
section_start = time.time()

# adata.write(results_file)
# adata

# PCA
sc.tl.pca(adata, svd_solver="arpack", n_comps=30)

print(f"[PROFILE] pca: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

# neighborhood graph
sc.pp.neighbors(adata, n_pcs=30)

print(f"[PROFILE] neighbors: {time.time() - section_start:.6f} s")


# %% 
# for fixing disconnected clusters or connectivity issues:
#sc.tl.paga(adata)
#sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#cs.tl.umap(adata, init_pos='paga')


# adata.write(results_file)
# adata

# %%
section_start = time.time()

# clustering
sc.tl.louvain(adata, resolution=0.5)

print(f"[PROFILE] louvain_clustering: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

# UMAP
sc.tl.umap(adata, n_components=30)

print(f"[PROFILE] umap: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

# write results
adata.write(results_file)

print(f"[PROFILE] write_h5ad: {time.time() - section_start:.6f} s")


# %%
section_start = time.time()

# Fine-grain profiling of the final rank_genes_groups step using cprofile.run() function with ranking function inside it
def run_rank_genes_groups():
    sc.tl.rank_genes_groups(
        adata,
        "louvain",
        method="wilcoxon",
        use_raw=True,
    )

cProfile.run("run_rank_genes_groups()", profile_file)

rank_elapsed = time.time() - section_start
print(f"[PROFILE] rank_genes_groups: {rank_elapsed:.6f} s")
print(f"[CPROFILE] saved: {profile_file}")
print(f"[OUTPUT] AnnData saved: {results_file}")
