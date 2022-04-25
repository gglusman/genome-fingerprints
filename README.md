# genome-fingerprints
Software for creating and comparing genome fingerprints.  
More information and datasets: http://db.systemsbiology.net/gestalt/genome_fingerprints/  
If you find Genome Fingerprints useful for your work, please cite:  
Glusman G, Mauldin DE, Hood L and Robinson M. Ultrafast comparison of personal genomes via precomputed genome fingerprints. Front. Genet. 2017 8:136.

To install this software, simply clone this repository:  
	`git clone https://github.com/gglusman/genome-fingerprints`

Here are example commands on how to compute fingerprints, and then how to compare them. Anything in **bold**, please replace with parameter values and file names that make sense for your project.

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

To compute fingerprints for a large dataset, like the Thousand Genomes Project (TGP), that is made available as a collection of per-chromosome multi-sample VCFs, use the following protocol. In the commands below, **TGP** represents the TGP version, **TGPdata** is the directory where you have your copy of the TGP data, and **TGPfp** represents the directory where the fingerprints will be stored.

1. Compute raw fingerprints per chromosome:  
	`bin/computeDMF1000genomes.pl` **TGPdata** **TGPfp**

2. Combine the raw fingerprints:  
	`bin/combineDMF1000genomes.pl` **TGPfp**

3. Normalize them:  
	`bin/normalizeDMFsInDir.pl` **TGPfp**/autosomal

4. Serialize the normalized fingerprints, using L=200:  
	`bin/serializeDMFs.pl` **TGPfp**/**TGP**-200 200 **TGPfp**/autosomal.norm/*.gz

5. Do all-against-all comparison within the set:  
	`bin/searchDMFs.pl` **TGPfp**/**TGP**-200 | sort -k3rn | gzip -c > **TGPfp**/**TGP**.aaa.gz

For example, if you have the GRCh37 version of the TGP in a directory called TGP37, the actual commands would be:  
`	bin/computeDMF1000genomes.pl TGP37 TGP37fp`  
`	bin/combineDMF1000genomes.pl TGP37fp`  
`	bin/normalizeDMFsInDir.pl TGP37fp/autosomal`  
`	bin/serializeDMFs.pl TGP37fp/TGP37-200 200 TGP37fp/autosomal.norm/*.gz`  
`	bin/searchDMFs.pl TGP37fp/TGP37-200 | sort -k3rn | gzip -c > TGP37fp/TGP37.aaa.gz`  
`	zcat TGP37fp/TGP37.aaa.gz | bin/findSurprises.pl`

Note that `bin/computeDMF1000genomes.pl` can be run in parallel copies to process more than one chromosome at a time, if your system has enough memory. Similarly, `bin/normalizeDMFsInDir.pl` can be run in concurrent jobs, if needed.
