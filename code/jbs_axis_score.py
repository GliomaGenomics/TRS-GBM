"""
Calculate JBS axis scores for an expression matrix.

Each sample's score is its log2(FPKM + pseudocount) profile dotted with the JBS
axis loadings; a patient's TRS is their recurrent tumour's score minus their
primary's. The matrix is genes × samples (-T for samples × genes), with gene
identifiers in the first column matching the gene metadata: Ensembl IDs, or
symbols with --no-Ensembl. Needs only pandas and numpy.

Morgan Thomas, Aug 2026.

Example:
    pixi run jbs_axis_score --mtx counts.csv --genes data/gene_metadata/gene_metadata_loadings.csv --sep "," -T
"""

import argparse
import sys

import numpy as np
import pandas as pd

LOADINGS_COLUMN = 'jbs_axis_loading'
LENGTH_COLUMN = 'length'


def build_parser():
    """Define the command line interface."""
    parser = argparse.ArgumentParser(
        description='Score samples on the JBS axis from an expression matrix and PCA gene loadings.'
    )

    parser.add_argument(
        '--mtx',
        required=True,
        help='path/to/expression_matrix.csv',
        type=str
    )
    parser.add_argument(
        '--genes',
        required=True,
        help='path/to/loadings.csv',
        type=str
    )
    parser.add_argument(
        '-o', '--out',
        default='axis_scores.txt',
        help='path/to/output.txt for the tab-separated scores. Default: axis_scores.txt',
        type=str
    )
    parser.add_argument(
        '--FPKM',
        default=True,
        action=argparse.BooleanOptionalAction,
        help='Whether to FPKM transform the input matrix (if providing a raw count matrix). Use --no-FPKM if the matrix is already normalised. Default: True',
    )
    parser.add_argument(
        '--pseudocount',
        default=0.01,
        help='Pseudocount added before the log2 transform, to keep zero counts finite. Default: 0.01',
        type=float
    )
    parser.add_argument(
        '--sep',
        default='\t',
        help=r'Separator for the count matrix. Default: "\t" (tab)',
        type=str
    )
    parser.add_argument(
        '-T',
        action='store_true',
        default=False,
        help='Transpose the count matrix, if it is shape (samples, genes).',
    )
    parser.add_argument(
        '--Ensembl',
        default=True,
        action=argparse.BooleanOptionalAction,
        help='Whether to use Ensembl IDs for gene identifiers in the first matrix column. Use --no-Ensembl to match on gene symbols instead. Default: True.',
    )

    return parser


IDENTIFIER_HELP = (
    "Check the correct gene identifiers are used, setting --no-Ensembl if gene symbols are used. "
    "Check the input matrix is shape (genes, samples), if (sample, gene) use the -T flag."
)


def fpkm(counts, lengths):
    """
    FPKM normalise a count matrix.

    Parameters:
    counts (pd.DataFrame): raw counts, shape (genes, samples).
    lengths (pd.Series): gene lengths in bases, indexed like the rows of counts.

    Returns:
    pd.DataFrame: FPKM values, shape (genes, samples).
    """
    library_sizes = counts.sum(axis=0) / 1e6
    rpm = counts.div(library_sizes, axis=1)
    return rpm.div(lengths / 1e3, axis=0)


def collapse_duplicates(count_matrix, gene_metadata):
    """Sum counts across duplicated gene identifiers, averaging their metadata."""
    if count_matrix.index.has_duplicates:
        n_duplicates = count_matrix.index.duplicated().sum()
        count_matrix = count_matrix.groupby(level=0).sum()
        print(f"Summed counts across {n_duplicates} duplicate gene identifiers in the count matrix.")

    if gene_metadata.index.has_duplicates:
        n_duplicates = gene_metadata.index.duplicated().sum()
        gene_metadata = gene_metadata.groupby(level=0).mean()
        print(f"Averaged length and loadings across {n_duplicates} duplicate gene identifiers in the gene metadata.")

    return count_matrix, gene_metadata


def axis_scores(count_matrix, gene_metadata, FPKM=True, pseudocount=0.01,
                length_column=LENGTH_COLUMN, loadings_column=LOADINGS_COLUMN):
    """
    Score every sample of a count matrix on the JBS axis.

    Parameters:
    count_matrix (pd.DataFrame): counts, shape (genes, samples).
    gene_metadata (pd.DataFrame): gene lengths and axis loadings, indexed by gene identifier.

    Returns:
    pd.Series: one axis score per sample.
    """
    gene_metadata = gene_metadata[[length_column, loadings_column]]
    n_loading_genes = gene_metadata[loadings_column].notna().sum()
    print(f"Loaded metadata for {gene_metadata.shape[0]} genes, "
          f"{n_loading_genes} of which carry an axis loading.")

    # collapse non-unique gene identifiers
    count_matrix, gene_metadata = collapse_duplicates(count_matrix, gene_metadata)

    # Restrict both to the shared gene set.
    n_genes = count_matrix.shape[0]
    intersection_genes = count_matrix.index.intersection(gene_metadata.index)

    if len(intersection_genes) == 0:
        sys.exit(
            "Error: no genes shared between the count matrix and the gene metadata.\n"
            + IDENTIFIER_HELP
        )

    count_matrix = count_matrix.loc[intersection_genes]
    gene_metadata = gene_metadata.loc[intersection_genes]
    print(f"Filtered count matrix from {n_genes} to {count_matrix.shape[0]} "
          f"({count_matrix.shape[0]/n_genes*100:.1f}%) genes present in the gene metadata.")

    if FPKM:
        # FPKM needs a length for every gene, so the library size is taken over
        # the genes that have one.
        missing_lengths = gene_metadata[length_column].isna()
        if missing_lengths.any():
            count_matrix = count_matrix.loc[~missing_lengths]
            gene_metadata = gene_metadata.loc[~missing_lengths]
            print(f"Dropped {missing_lengths.sum()} genes with no recorded length.")
        expression = fpkm(count_matrix, gene_metadata[length_column])
    else:
        expression = count_matrix
        print("Skipping FPKM transform, assuming the input matrix is already normalised.")

    # Log2 + pseudocount
    log_expression = np.log2(expression + pseudocount)

    # Dot the log expression (genes, samples) with the loadings (genes,) to give a score per sample.
    # The PCA was fit with center = FALSE and scale = FALSE, so no standardisation is applied here.
    loadings = gene_metadata[loadings_column].dropna()
    print(f"Retained {loadings.shape[0]} of {n_loading_genes} "
          f"({loadings.shape[0]/n_loading_genes*100:.1f}%) loading genes after matching and filtering.")
    scores = log_expression.loc[loadings.index].T.dot(loadings)
    print(f"Scored {scores.shape[0]} samples using {loadings.shape[0]} genes with axis loadings.")
    return scores


def main(args):
    """Score every sample in the input matrix and write the results to args.out."""
    # Load count/FPKM matrix
    count_matrix = pd.read_csv(args.mtx, sep=args.sep, index_col=0)
    if args.T:
        count_matrix = count_matrix.T
    # A gene left unquantified in a sample counts as unexpressed.
    count_matrix = count_matrix.fillna(0)

    # Load PCA loadings and gene lengths for FPKM
    gene_metadata = pd.read_csv(args.genes, index_col=0)
    if not args.Ensembl:
        gene_metadata = gene_metadata.set_index('GeneName')

    scores = axis_scores(count_matrix, gene_metadata,
                         FPKM=args.FPKM, pseudocount=args.pseudocount)

    scores.to_frame('axis_score').rename_axis('Sample').to_csv(args.out, sep='\t')
    print(f"Wrote scores to {args.out}")


if __name__ == '__main__':
    main(build_parser().parse_args())
