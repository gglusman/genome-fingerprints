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

This project is related to (but distinct from) the Genotype Fingerprints: https://github.com/gglusman/genotype-fingerprints
