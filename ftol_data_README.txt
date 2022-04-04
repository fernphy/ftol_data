This ftol_data_README.txt file was generated on 2022-04-04 by Joel Nitta

--------------------------------------------------------------------------------

GENERAL INFORMATION

--------------------------------------------------------------------------------

Title of Dataset: Fern Tree of Life (FTOL) and associated data

Principal Investigator: Joel H. Nitta

Department of Biological Sciences, Graduate School of Science, The University of
Tokyo, 2-11-16 Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan joelnitta@gmail.com

Associate or Co-investigators: Eric Schuettpelz, Santiago Ramírez-Barahona,
Wataru Iwasaki

Date of data collection: 1990 - 2021

Geographic location of data collection: Global

Information about funding sources that supported the collection of the data:
Funding provided in part by the Japan Society for the Promotion of Science
(Kakenhi) Grant Number 16H06279 and the Smithsonian National Museum of Natural
History Peter Buck Fellowship (JHN).

--------------------------------------------------------------------------------

SHARING/ACCESS INFORMATION

--------------------------------------------------------------------------------

Licenses/restrictions placed on the data: CC0 v1.0 license

Links to publications that cite or use the data:

Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. 2022. An open and
continuously updated fern tree of life (FTOL).
https://doi.org/10.1101/2022.03.31.486640

Links to other publicly accessible locations of the data: none.

Links/relationships to ancillary data sets:

-   ferncal (https://github.com/fernphy/ferncal)
-   pteridocat (https://github.com/fernphy/pteridocat)
-   FTOL input data (https://doi.org/10.6084/m9.figshare.19474316.v1)

Was data derived from another source? Yes, in part from GenBank
(https://www.ncbi.nlm.nih.gov/genbank/), which places no restrictions on its use
or distribution.

Recommended citation for this dataset:

FTOL Working Group (2022). Fern Tree of Life (FTOL) input data. FIXME: add DOI
when available

--------------------------------------------------------------------------------

DATA & FILE OVERVIEW

--------------------------------------------------------------------------------

File List:

-   ftol_acc_table_long.csv: GenBank accessions used in the FTOL, long format.
-   ftol_acc_table_wide.csv: GenBank accessions used in FTOL, wide format.
-   ftol_plastome_alignment.fasta.gz: Aligned plastome DNA sequences used to
    build FTOL.
-   ftol_plastome_con.tre: FTOL backbone phylogeny.
-   ftol_plastome_parts.csv: Start and end positions of loci in plastome DNA
    sequence alignment.
-   ftol_sanger_alignment.fasta.gz: Aligned (mostly) Sanger DNA sequences used
    to build FTOL.
-   ftol_sanger_con_dated.tre: FTOL dated consensus phylogeny.
-   ftol_sanger_con_fossils.csv: Fossil calibration points used for dating FTOL
    consensus phylogeny.
-   ftol_sanger_con.tre: FTOL consensus phylogeny.
-   ftol_sanger_ml_dated.tre: FTOL dated maximum-likelihood phylogeny.
-   ftol_sanger_ml_fossils.csv: Fossil calibration points used for dating FTOL
    ML phylogeny.
-   ftol_sanger_ml.tre: FTOL maximum-likelihood phylogeny.
-   ftol_sanger_parts.csv: Start and end positions of loci in Sanger DNA
    sequence alignment.
-   ftol_sanger_sampling.csv: Taxonomic data of species in FTOL.

--------------------------------------------------------------------------------

METHODOLOGICAL INFORMATION

--------------------------------------------------------------------------------

The data included are the results of a pipeline that (mostly) automatically
generates a maximally sampled fern phylogenetic tree based on plastid sequences
in GenBank (https://github.com/fernphy/ftol).

The first step is to generate a set of reference FASTA files for 79 target loci
(one per locus). These include 77 protein-coding genes based on a list of 83
genes (Wei et al. 2017) that was filtered to only genes that show no evidence of
duplication, plus two spacer regions (trnL-trnF and rps4-trnS). This is done
with custom R scripts contained in https://github.com/fernphy/ftol, in
particular prep_ref_seqs_plan.R
(https://github.com/fernphy/ftol/blob/main/prep_ref_seqs_plan.R).

Next, all available fern accessions for seven target “Sanger loci” (plastid
regions typically sequenced using Sanger technology) and all available fern
plastomes (accessions >7000 bp) were downloaded from GenBank. Selected non-fern
accessions (outgroups) were downloaded as well. Sequences matching the target
loci were then extracted from each accesion using the reference FASTA files with
the “Reference_Blast_Extract.py” script of superCRUNCH (Portik and Wiens 2020).
Putative rogues (i.e., misidentifications or contaminations) were identified and
removed after sequence extraction either by manual inspection or all-by-all
BLAST (Altschul et al. 1997).

The extracted sequences were aligned with MAFFT (Katoh et al. 2002).
Phylogenetic analysis was done in two steps with maximum-likelihood in IQ-TREE
(Nguyen et al. 2015). First, a backbone phylogeny (ftol_plastome_con.tre) was
inferred using a matrix from whole plastome sequences
(ftol_plastome_alignment.fasta.gz). Next, the backbone tree was used as a
constraint tree in ML analysis of the “Sanger” dataset (loci commonly obtained
by Sanger sequencing; ftol_sanger_alignment.fasta.gz); this resulted in a ML
tree (ftol_sanger_ml.tre) and the extended majority-rule consenus of 1,000
bootstrap trees (ftol_sanger_con.tre). Molecular dating analysis was carried out
separately on the ML tree and consensus tree each using a set of fossil
calibration points (ftol_sanger_ml_fossils.csv and ftol_sanger_con_fossils.csv,
respectively), resulting in ultrametric trees (ftol_sanger_ml_dated.tre and
ftol_sanger_con_dated.tre, respectively).

For additional methodological details, see:

Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. 2022. An open and
continuously updated fern tree of life (FTOL).
https://doi.org/10.1101/2022.03.31.486640

--------------------------------------------------------------------------------

DATA-SPECIFIC INFORMATION

--------------------------------------------------------------------------------

ftol_acc_table_long.csv: GenBank accessions used in the FTOL, long format.

Number of variables: 8

Number of cases/rows: 14651

Variable list:

-   species: Species name; matches names of tips in tree
-   locus: Name of locus (gene or intergenic spacer region)
-   accession: GenBank accession number
-   seq_len: Sequence length (bp), excluding any missing or ambiguous bases
-   sci_name: Scientific name used in FTOL
-   ncbi_name: Scientific name used in the NCBI taxonomic database
-   ncbi_taxid: NCBI taxonomy database unique identifier
-   outgroup: Logical; TRUE for outgroup taxa, FALSE for ingroup taxa (ferns)

Missing data codes: No missing data.

Specialized formats or other abbreviations used: None.

MD5 checksum: 76047ae746a7a1bbf7cb11316ca13867

--------------------------------------------------------------------------------

ftol_acc_table_wide.csv: GenBank accessions used in FTOL, wide format.

Number of variables: 13

Number of cases/rows: 5583

Variable list:

-   species: Species name; matches names of tips in tree
-   atpA: GenBank accession number for atpA
-   atpB: GenBank accession number for atpB
-   matK: GenBank accession number for matK
-   rbcL: GenBank accession number for rbcL
-   rps4: GenBank accession number for rps4
-   rps4-trnS: GenBank accession number for rps4-trnS
-   trnL-trnF: GenBank accession number for trnL-trnF
-   plastome: GenBank accession number for plastomes
-   join_by: Method used to join loci
-   specimen_voucher: Specimen voucher
-   publication: Publication
-   outgroup: Logical; TRUE for outgroup taxa, FALSE for ingroup taxa (ferns)

Missing data codes: ‘NA’ for missing or inapplicable data.

Specialized formats or other abbreviations used: None.

MD5 checksum: 0e818a995aeaaed5b411eb7045e4616e

--------------------------------------------------------------------------------

ftol_plastome_alignment.fasta.gz: Aligned plastome DNA sequences used to build
FTOL. In compressed (tar.gz) FASTA format. Includes 79 concatenated loci. The
start and end position (column) of each locus is given in
ftol_plastome_parts.csv.

Number of bases (columns): 74674

Number of rows (taxa): 429

Specialized formats or other abbreviations used: None.

MD5 checksum: da6946973c5e2a4334892790206e0a4f

--------------------------------------------------------------------------------

ftol_plastome_con.tre: FTOL backbone phylogeny. Inferred using
maximum-likelihood from DNA sequences in ftol_plastome_alignment.fasta.gz.
Extended majority-rule consensus of 1000 bootstrap trees. Rooted on bryophytes.
In newick format.

Number of tips: 429

Specialized formats or other abbreviations used: None.

MD5 checksum: 24a3293bb0e8cf5bf095558b951561a7

--------------------------------------------------------------------------------

ftol_plastome_parts.csv: Start and end positions of loci in plastome DNA
sequence alignment.

Number of variables: 3

Number of cases/rows: 79

Variable list:

-   locus: Name of locus (gene or intergenic spacer region)
-   start: Start position (column number) of locus in concatenated DNA alignment
-   end: End position (column number) of locus in concatenated DNA alignment

Missing data codes: None.

Specialized formats or other abbreviations used: None.

MD5 checksum: c0a11e1965de8e4c46f697da1c2607a8

--------------------------------------------------------------------------------

ftol_sanger_alignment.fasta.gz: Aligned (mostly) Sanger DNA sequences used to
build FTOL. In compressed (tar.gz) FASTA format. Includes 7 concatenated loci.
The start and end position (column) of each locus is given in
ftol_sanger_parts.csv.

Number of bases (columns): 12118

Number of rows (taxa): 5581

Specialized formats or other abbreviations used: None.

MD5 checksum: 61a60db4d2d412e2f76af2a5a157836c

--------------------------------------------------------------------------------

ftol_sanger_con_dated.tre: FTOL dated consensus phylogeny. Inferred using
maximum-likelihood from DNA sequences in ftol_plastome_alignment.fasta.gz.
Extended majority-rule consensus of 1000 bootstrap trees. Rooted on bryophytes.
Divergence times estimated with fossil calibration points
(ftol_sanger_con_fossils.csv) using treePL. In newick format.

Number of tips: 5581

Specialized formats or other abbreviations used: None.

MD5 checksum: baf2d8b00c46ece96e1af350ac56c4e0

--------------------------------------------------------------------------------

ftol_sanger_con_fossils.csv: Fossil calibration points used for dating FTOL
consensus phylogeny. A subset of fossil data contained in the ‘ferncal’ fossil
database, v.1.0.0 (https://doi.org/10.5281/zenodo.6395323).

The fossils in ‘ftol_sanger_con_fossils.csv’ and ‘ftol_sanger_ml_fossils.csv’
are the same, but the node each calibrates may differ between the trees because
of differences in topology.

The node corresponding to the fossil constraint is defined as the most recent
common ancestor (MRCA, column ‘mrca’) of two tips columns (‘tip_1’ and ‘tip_2’)
for crown affinities, or its parent node (column ‘stem_mrca’) for stem
affinities. The two tips are identified automatically for monophyletic clades,
or by hand for non-monophyletic clades. ‘mrca’ is not defined for monotypic
groups (only ‘stem_mrca’).

Does not include the constraint on the root of the tree (landplants; 475 Ma).

Number of variables: 12

Number of cases/rows: 66

Variable list:

-   n_fos: Unique ID number for fossil
-   minimum_age: Minimum age to apply to fossil constraint
-   node_calibrated: Node calibrated by fossil constraint. Combination of
    ‘affinities’ and ‘affinities_group’
-   fossil_taxon: Taxonomic name of fossil (without author)
-   affinities_group: Type of group the fossil belongs to (crown or stem)
-   affinities: Narrowest clade the fossil belongs to; the clade whose date is
    constrained by the fossil
-   monophyly: Are the affinities monophyletic? ‘Yes,’ ‘No,’ or ‘Monotypic’
-   number_tips: Number of tips in the clade constrained by the fossil
-   mrca: Node number of MRCA for the clade constrained by the fossil
-   stem_mrca: Node number of the parent node of the MRCA for the clade
    constrained by the fossil
-   tip_1: Name of one taxon that defines the clade constrained by the fossil
-   tip_2: Name of another taxon that defines the clade constrained by the
    fossil

Missing data codes: ‘NA’ for missing or inapplicable data.

Specialized formats or other abbreviations used: None.

MD5 checksum: bab06778b3cb320cb2664c0515205eb0

--------------------------------------------------------------------------------

ftol_sanger_con.tre: FTOL consensus phylogeny. Inferred using maximum-likelihood
from DNA sequences in ftol_plastome_alignment.fasta.gz. Extended majority-rule
consensus of 1000 bootstrap trees. Rooted on bryophytes. In newick format.

Number of tips: 5581

Specialized formats or other abbreviations used: None.

MD5 checksum: 09890cc6b8c752782b306324031f6125

--------------------------------------------------------------------------------

ftol_sanger_ml_dated.tre: FTOL dated maximum-likelihood phylogeny. Inferred
using maximum-likelihood from DNA sequences in ftol_plastome_alignment.fasta.gz.
Rooted on bryophytes. Divergence times estimated with fossil calibration points
(ftol_sanger_con_fossils.csv) using treePL. In newick format.

Number of tips: 5581

Specialized formats or other abbreviations used: None.

MD5 checksum: e806b8d346e7fb0c77b06c32fd6b8fca

--------------------------------------------------------------------------------

ftol_sanger_ml_fossils.csv: Fossil calibration points used for dating FTOL
maximum-likelihood phylogeny. A subset of fossil data contained in the ‘ferncal’
fossil database, v.1.0.0 (https://doi.org/10.5281/zenodo.6395323).

The fossils in ‘ftol_sanger_con_fossils.csv’ and ‘ftol_sanger_ml_fossils.csv’
are the same, but the node each calibrates may differ between the trees because
of differences in topology.

For more details, see entry for ftol_sanger_con_fossils.csv.

Does not include the constraint on the root of the tree (landplants; 475 Ma).

Number of variables: 12

Number of cases/rows: 66

Variable list: See entry for ftol_sanger_con_fossils.csv

Missing data codes: ‘NA’ for missing or inapplicable data.

Specialized formats or other abbreviations used: None.

MD5 checksum: bab06778b3cb320cb2664c0515205eb0

--------------------------------------------------------------------------------

ftol_sanger_ml.tre: FTOL maximum-likelihood phylogeny. Inferred using
maximum-likelihood from DNA sequences in ftol_plastome_alignment.fasta.gz.
Rooted on bryophytes. In newick format.

Number of tips: 5581

Specialized formats or other abbreviations used: None.

MD5 checksum: 09890cc6b8c752782b306324031f6125

--------------------------------------------------------------------------------

ftol_sanger_parts.csv: Start and end positions of loci in plastome DNA sequence
alignment.

Number of variables: 3

Number of cases/rows: 7

Variable list:

-   locus: Name of locus (gene or intergenic spacer region)
-   start: Start position (column number) of locus in concatenated DNA alignment
-   end: End position (column number) of locus in concatenated DNA alignment

Missing data codes: None.

Specialized formats or other abbreviations used: None.

MD5 checksum: e70315c2c36169da38cbac7024a6c0d7

--------------------------------------------------------------------------------

ftol_sanger_sampling.csv: Taxonomic data of species in FTOL.

Number of variables: 9

Number of cases/rows: 5581

Variable list:

-   species: Species name
-   genus: Genus name
-   order: Order name
-   suborder: Suborder name
-   family: Family name
-   subfamily: Subfamily name
-   major_clade: Informal higher level clade name, either order or suborder
-   outgroup: Logical; TRUE for outgroup taxa, FALSE for ingroup taxa (ferns)
-   subgenus: Subgenus name (only provided if used for molecular dating)

Missing data codes: ‘NA’ for missing or inapplicable data.

Specialized formats or other abbreviations used: None.

MD5 checksum: dee9b1b2f3cc5322a775470a0d1b5832

--------------------------------------------------------------------------------

CHANGE LOG

2022-04-04

-   Add DOI for Nitta et al. 2022 preprint
-   Change name of README file from README.txt to ftol_data_README.txt

2022-03-31

-   Generate this README file.

--------------------------------------------------------------------------------

REFERENCES

Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang,
Zheng Zhang, Webb Miller, and David J. Lipman. 1997. “Gapped BLAST and
PSI-BLAST: A New Generation of Protein Database Search Programs.” Nucleic Acids
Research 25: 3389–3402. https://doi.org/10.1093/nar/25.17.3389.

Katoh, Kazutaka, Kazuharu Misawa, Keiichi Kuma, and Takashi Miyata. 2002.
“MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based on Fast
Fourier Transform.” Nucleic Acids Research 30 (14): 3059–66.
https://doi.org/10.1093/nar/gkf436.

Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh.
2015. “IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating
Maximum-Likelihood Phylogenies.” Molecular Biology and Evolution 32 (1): 268–74.
https://doi.org/10.1093/molbev/msu300.

Portik, Daniel M., and John J. Wiens. 2020. “SuperCRUNCH: A Bioinformatics
Toolkit for Creating and Manipulating Supermatrices and Other Large Phylogenetic
Datasets.” Edited by David Orme. Methods in Ecology and Evolution 11 (6):
763–72. https://doi.org/ggx588.

Wei, Ran, Yue-Hong Yan, AJ Harris, Jong-Soo Kang, Hui Shen, Qiao-Ping Xiang, and
Xian-Chun Zhang. 2017. “Plastid Phylogenomics Resolve Deep Relationships Among
Eupolypod II Ferns with Rapid Radiation and Rate Heterogeneity.” Genome Biology
and Evolution 9 (6): 1646–57. https://doi.org/10.1093/gbe/evx107.

Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang,
Zheng Zhang, Webb Miller, and David J. Lipman. 1997. “Gapped BLAST and
PSI-BLAST: A New Generation of Protein Database Search Programs.” Nucleic Acids
Research 25: 3389–3402. https://doi.org/10.1093/nar/25.17.3389.

Katoh, Kazutaka, Kazuharu Misawa, Keiichi Kuma, and Takashi Miyata. 2002.
“MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based on Fast
Fourier Transform.” Nucleic Acids Research 30 (14): 3059–66.
https://doi.org/10.1093/nar/gkf436.

Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh.
2015. “IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating
Maximum-Likelihood Phylogenies.” Molecular Biology and Evolution 32 (1): 268–74.
https://doi.org/10.1093/molbev/msu300.

Portik, Daniel M., and John J. Wiens. 2020. “SuperCRUNCH: A Bioinformatics
Toolkit for Creating and Manipulating Supermatrices and Other Large Phylogenetic
Datasets.” Edited by David Orme. Methods in Ecology and Evolution 11 (6):
763–72. https://doi.org/ggx588.

Wei, Ran, Yue-Hong Yan, AJ Harris, Jong-Soo Kang, Hui Shen, Qiao-Ping Xiang, and
Xian-Chun Zhang. 2017. “Plastid Phylogenomics Resolve Deep Relationships Among
Eupolypod II Ferns with Rapid Radiation and Rate Heterogeneity.” Genome Biology
and Evolution 9 (6): 1646–57. https://doi.org/10.1093/gbe/evx107.
