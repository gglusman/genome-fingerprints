# genome-fingerprints
Software for creating and comparing genome fingerprints.  
More information and datasets: http://db.systemsbiology.net/gestalt/genome_fingerprints/  
If you find Genome Fingerprints useful for your work, please cite:  
Glusman G, Mauldin DE, Hood L and Robinson M. Ultrafast comparison of personal genomes via precomputed genome fingerprints. Front. Genet. 2017 8:136.

## Installing

To install this software, simply clone this repository:  
	`git clone https://github.com/gglusman/genome-fingerprints`

This software works best on whole-genome sequencing (WGS) data, in VCF and some other formats. The data type analyzed by this software is variants relative to a reference genome. Therefore, FASTQ, BAM, CRAM and the like are not valid input formats.

## Computing and comparing fingerprints

Here are example commands on how to compute fingerprints from genomes that are represented separately, each one in its own VCF, and then how to compare the resulting fingerprints. Anything in **bold**, please replace with parameter values and file names that make sense for your project.

1. To create a fingerprint for a genome:  
	`bin/computeDMF.pl` **myGenome** **path-to-my-vcfs/myGenome.vcf.gz**  
	...will generate myGenome.outn.gz (and some other files)

2. To compare two fingerprints:  
	`bin/compareDMFs.pl` **myFirstGenome.outn.gz** **mySecondGenome.outn.gz**

3. To serialize fingerprints into a database, using L=120:  
	`bin/serializeDMFs.pl` **myFingerprintCollection** 120 @**myListOfFingerprints**  
	`bin/serializeDMFs.pl` **myFingerprintCollection** 120 **\*.outn.gz**

4. To compare a fingerprint to a database:  
	`bin/searchDMFs.pl` **myGenome.outn.gz** **myFingerprintCollection**  
	...see the data directory for an example database (CEPH1463 pedigree)

5. To compare two databases:  
	`bin/searchDMFs.pl` **aFingerprintCollection** **anotherFingerprintCollection** > **setAvsB.aaa**

6. To perform all-against-all comparisons in one database:  
	`bin/searchDMFs.pl` **aFingerprintCollection** > **setA.aaa**

7. To find surprises in the comparison:  
	cat **comparison.aaa** | `bin/findSurprises.pl` > **surprises.txt** 

## Working with large datasets

To compute fingerprints for a large dataset, like the Thousand Genomes Project (TGP), that is made available as a collection of per-chromosome multi-sample VCFs, use the following protocol. In the commands below, **TGP** represents the TGP version, **TGPdata** is the directory where you have your copy of the TGP data, and **TGPfp** represents the directory where the fingerprints will be stored.

1. Compute raw fingerprints per chromosome, using L=120 and L=200:  
	`bin/computeDMF1000genomes.pl` **TGPdata** 120,200 **TGPfp**

2. Combine the raw fingerprints:  
	`bin/combineDMF1000genomes.pl` **TGPfp**

3. Normalize them:  
	`bin/normalizeDMFsInDir.pl` **TGPfp**/autosomal

4. Serialize the normalized fingerprints, using L=200:  
	`bin/serializeDMFs.pl` **TGPfp**/**TGP**-200 200 **TGPfp**/autosomal.norm/*.gz

5. Do all-against-all comparison within the set:  
	`bin/searchDMFs.pl` **TGPfp**/**TGP**-200 | sort -k3rn | gzip -c > **TGPfp**/**TGP**.aaa.gz

For example, if you have the GRCh37 version of the TGP in a directory called TGP37, and you wanted to produce and compare fingerprints with L=200, the actual commands would be:  
`	bin/computeDMF1000genomes.pl TGP37 200 TGP37fp`  
`	bin/combineDMF1000genomes.pl TGP37fp`  
`	bin/normalizeDMFsInDir.pl TGP37fp/autosomal`  
`	bin/serializeDMFs.pl TGP37fp/TGP37-200 200 TGP37fp/autosomal.norm/*.gz`  
`	bin/searchDMFs.pl TGP37fp/TGP37-200 | sort -k3rn | gzip -c > TGP37fp/TGP37.aaa.gz`  
`	zcat TGP37fp/TGP37.aaa.gz | bin/findSurprises.pl`

Note that `bin/computeDMF1000genomes.pl` can be run in parallel copies to process more than one chromosome at a time, if your system has enough memory. Similarly, `bin/normalizeDMFsInDir.pl` can be run in concurrent jobs, if needed.

The TGP datasets can be downloaded as follows. Beware, these are big downloads!

TGP37: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/  
TGP38L: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/  
TGP38C: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/  
TGP38S: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/  
TGP38X: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/  
TGP38H: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/  
TGP38N and TGP38Nr: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/  
TGP37r: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/related_samples_vcf/  
TGP38Sr: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/supporting/related_samples/  
TGP38Xr: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/supporting/related_samples/

