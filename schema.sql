
CREATE TABLE Annotation (
	Key     	   	  VARCHAR(50) NOT NULL,
	File     		  VARCHAR(50) NOT NULL,
    Chrom    		  VARCHAR(50) NOT NULL,
    Start    	      VARCHAR(50) NOT NULL,
    REF      		  VARCHAR(50) NOT NULL,
    ALT     		  VARCHAR(50) NOT NULL,
	Type    		  VARCHAR(10) NOT NULL,
    Gene     		  VARCHAR(50) NOT NULL,
    Feature  	      VARCHAR(50) NOT NULL,
    FeatureType       VARCHAR(50) NOT NULL,
    Consequence       VARCHAR(100),
    cDNAPosition      VARCHAR(100),
    CDSPosition       VARCHAR(100),
    ProteinPosition   VARCHAR(100),
    AminoAcids 	      VARCHAR(100),
    Codons  		  VARCHAR(100),
    ExistingVariation VARCHAR(100),
	Extra			  VARCHAR(500)
);





/*

CREATE TABLE Details (
	TransName     VARCHAR(50) NOT NULL,
	SYMBOL        VARCHAR(50), -- Gene symbol (e.g. HGNC) 
	SYMBOL_SOURCE VARCHAR(50), -- Source of gene symbol
	HGNC_ID       VARCHAR(50), -- Stable identifer of HGNC gene symbol
	BIOTYPE       VARCHAR(50), -- Biotype of transcript or regulatory feature
	CANONICAL     VARCHAR(50), -- Indicates if transcript is canonical for this gene
	TSL	      	  VARCHAR(50), -- Transcript support level
	APPRIS	      VARCHAR(50), -- Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods
	CCDS	      VARCHAR(50), -- Indicates if transcript is a CCDS transcript
	ENSP	      VARCHAR(50), -- Protein identifer
	SWISSPROT	  VARCHAR(50), -- UniProtKB/Swiss-Prot accession
	TREMBL	 	  VARCHAR(50), -- UniProtKB/TrEMBL accession
	UNIPARC	  	  VARCHAR(50), -- UniParc accession
	GENE_PHENO	  VARCHAR(50), -- Indicates if gene is associated with a phenotype, disease or trait
	SIFT	 	  VARCHAR(50), -- SIFT prediction and/or score
	PolyPhen	  VARCHAR(50), -- PolyPhen prediction and/or score
	EXON	 	  VARCHAR(50), -- Exon number(s) / total
	INTRON	 	  VARCHAR(50), -- Intron number(s) / total
	DOMAINS	  	  VARCHAR(50), -- The source and identifer of any overlapping protein domains
	miRNA	 	  VARCHAR(50), -- SO terms of overlapped miRNA secondary structure feature(s)
	AF	 	      VARCHAR(50), -- Frequency of existing variant in 1000 Genomes combined population
	AFR_AF	 	  VARCHAR(50), -- Frequency of existing variant in 1000 Genomes combined African population
	AMR_AF	 	  VARCHAR(50), -- Frequency of existing variant in 1000 Genomes combined American population
	EAS_AF	 	  VARCHAR(50), -- Frequency of existing variant in 1000 Genomes combined East Asian population
	EUR_AF	 	  VARCHAR(50), -- Frequency of existing variant in 1000 Genomes combined European population
	SAS_AF	 	  VARCHAR(50), -- Frequency of existing variant in 1000 Genomes combined South Asian population
	AA_AF	 	  VARCHAR(50), -- Frequency of existing variant in NHLBI-ESP African American population
	EA_AF	 	  VARCHAR(50), -- Frequency of existing variant in NHLBI-ESP European American population	
	gnomAD_AF	  VARCHAR(50), -- Frequency of existing variant in gnomAD exomes combined population
	gnomAD_AFR_AF VARCHAR(50), -- Frequency of existing variant in gnomAD exomes African/American population
	gnomAD_AMR_AF VARCHAR(50), -- Frequency of existing variant in gnomAD exomes American population
	gnomAD_ASJ_AF VARCHAR(50), -- Frequency of existing variant in gnomAD exomes Ashkenazi Jewish population
	gnomAD_EAS_AF VARCHAR(50), -- Frequency of existing variant in gnomAD exomes East Asian population
	gnomAD_FIN_AF VARCHAR(50), -- Frequency of existing variant in gnomAD exomes Finnish population
	gnomAD_NFE_AF VARCHAR(50), -- Frequency of existing variant in gnomAD exomes Non-Finnish European population
	gnomAD_OTH_AF VARCHAR(50), -- Frequency of existing variant in gnomAD exomes other combined populations
	gnomAD_SAS_AF VARCHAR(50), -- Frequency of existing variant in gnomAD exomes South Asian population	
	MAX_AF	 	  VARCHAR(50), -- Maximum observed allele frequency in 1000 Genomes, ESP and ExAC/gnomAD
	MAX_AF_POPS	  VARCHAR(50), -- Populations in which maximum allele frequency was observed
	CLIN_SIG	  VARCHAR(50), -- ClinVar clinical significance of the dbSNP variant
	SOMATIC	 	  VARCHAR(50), -- Somatic status of existing variant
	PHENO	 	  VARCHAR(50), -- Indicates if existing variant(s) is associated with a phenotype, disease or trait; multiple values correspond to multiple variants
	PUBMED	 	  VARCHAR(50), -- Pubmed ID(s) of publications that cite existing variant
	MOTIF_NAME	  VARCHAR(50), -- The source and identifier of a transcription factor binding profile (TFBP) aligned at this position
	MOTIF_POS	  VARCHAR(50), -- The relative position of the variation in the aligned TFBP
	HIGH_INF_PO   VARCHAR(50), -- A flag indicating if the variant falls in a high information position of the TFBP
	MOTIF_SCORE_CHANGE VARCHAR(50), -- The difference in motif score of the reference and variant sequences for the TFBP
	
    FOREIGN KEY (TransName)  REFERENCES Transcript(Name)
);
*/

