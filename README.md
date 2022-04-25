# genome-fingerprints
Software for creating and comparing genome fingerprints.  
More information and datasets: http://db.systemsbiology.net/gestalt/genome_fingerprints/  
If you find Genome Fingerprints useful for your work, please cite:  
Glusman G, Mauldin DE, Hood L and Robinson M. Ultrafast comparison of personal genomes via precomputed genome fingerprints. Front. Genet. 2017 8:136.

1. To create a fingerprint for a genome:  
	`bin/computeDMF.pl` _myGenome path-to-my-vcfs/myGenome.vcf.gz_  
	...will generate myGenome.outn.gz (and some other files)

2. To compare two fingerprints:  
	`bin/compareDMFs.pl` _myFirstGenome.outn.gz mySecondGenome.outn.gz_

3. To serialize fingerprints into a database, using L=120:  
	`bin/serializeDMFs.pl` _myFingerprintCollection_ 120 _@myListOfFingerprints_  
	`bin/serializeDMFs.pl` _myFingerprintCollection_ 120 _*.outn.gz_

4. To compare a fingerprint to a database:  
	`bin/searchDMFs.pl` _myGenome.outn.gz myFingerprintCollection_  
	...see the data directory for an example database (CEPH1463 pedigree)

5. To compare two databases:  
	`bin/searchDMFs.pl` _aFingerprintCollection anotherFingerprintCollection_

6. To perform all-against-all comparisons in one database:  
	`bin/searchDMFs.pl` _aFingerprintCollection_

To compute fingerprints for a large dataset, like the Thousand Genomes Project (TGP), that is made available as a collection of per-chromosome multi-sample VCFs, use the following protocol. In the commands below, _TGP_ represents the TGP version, _TGPdata_ is the directory where you have your copy of the TGP data, and _TGPfp_ represents the directory where the fingerprints will be stored.

1. Compute raw fingerprints per chromosome:  
	`bin/computeDMF1000genomes.pl` _TGPdata_ _TGPfp_

2. Combine the raw fingerprints:  
	`bin/combineDMF1000genomes.pl` _TGPfp_

3. Normalize them:  
	`bin/normalizeDMFsInDir.pl` _TGPfp_/autosomal

4. Serialize the normalized fingerprints, using L=200:  
	`bin/serializeDMFs.pl` _TGPfp_/_TGP_-200 200 _TGPfp_/autosomal.norm/*.gz

5. Do all-against-all comparison within the set:  
	`bin/searchDMFs.pl` _TGPfp_/_TGP_-200 | sort -k3rn | gzip -c > _TGPfp_/_TGP_.aaa.gz

For example, if you have the GRCh37 version of the TGP in a directory called TGP37, the actual commands would be:
`	bin/computeDMF1000genomes.pl TGP37 TGP37fp  
	bin/combineDMF1000genomes.pl TGP37fp  
	bin/normalizeDMFsInDir.pl TGP37fp/autosomal  
	bin/serializeDMFs.pl TGP37fp/TGP37-200 200 TGP37fp/autosomal.norm/*.gz  
	bin/searchDMFs.pl TGP37fp/TGP37-200 | sort -k3rn | gzip -c > TGP37fp/TGP37.aaa.gz`

Note that `bin/computeDMF1000genomes.pl` can be run in parallel copies to process more than one chromosome at a time, if your system has enough memory. Similarly, `bin/normalizeDMFsInDir.pl` can be run in concurrent jobs, if needed.
