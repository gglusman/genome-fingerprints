# genome-fingerprints
Software for creating and comparing genome fingerprints.
More information and datasets: http://db.systemsbiology.net/gestalt/genome_fingerprints/

To create a fingerprint for a genome:
	bin/computeDMF.pl myGenome path-to-my-vcfs/myGenome.vcf.gz
	...will generate myGenome.outn.gz (and some other files)

To compare two fingerprints:
	bin/compareDMFs.pl myFirstGenome.outn.gz mySecondGenome.outn.gz

To serialize fingerprints into a database, using L=120:
	bin/serializeDMFs.pl myFingerprintCollection 120 @myListOfFingerprints
	bin/serializeDMFs.pl myFingerprintCollection 120 *.outn.gz

To compare a fingerprint to a database:
	bin/searchDMFs.pl myFirstGenome.outn.gz myFingerprintCollection
	...see the data directory for an example database (CEPH1463 pedigree)

To compare two databases:
	bin/searchDMFs.pl aFingerprintCollection anotherFingerprintCollection

To perform all-against-all comparisons in one database:
	bin/searchDMFs.pl aFingerprintCollection

