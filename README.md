# Accelerating K-Means Clustering for Documents<br> with an Architecture-Friendly Pruning Method
This is a code set for appying our K-means clustering algorithm, AF-ICP, 
to **Large-scale and High-dimensional sparse data sets** like documents.
The codes are implemented with C.

## Requirements for executing codes
1. OS: CentOS 7 and later
2. g++ (GCC): >= 8.2.0
3. perl: >= 5.16
4. perf: 3.10

## Quick start: AF-ICP in five iterations
1. Prepare the 8.2M-sized PubMed data set with a procedure in [dataset](./dataset).<br>
This procedure creates ``./dataset/8.2M_pubmed.db``.
2. Execute ``make -f Makefile_itr5_aficp`` in ``./src``.<br>
 This makes ``./bin/itr5_aficp`` object in your system.
3. Execute the perl script ``./itr5_exeAFICP_8.2Mpubmed_perf.pl`` in ``./exe``.<br>
 The 8.2M-sized PubMed data set is loaded from ``./dataset/8.2M_pubmed.db`` (3.8GB)
 in around two minutes
 and given K=10,000, AF-ICP is executed with 50-thread parallel processing.<br>
 A log file is generated in ``./Log``.<br>

## Compare AF-ICP with other algorithms, ICP, TA-ICP, and CS-ICP
Go to [Comparison](./Comparison).

## License
Please check [LICENSE](LICENSE_v2.1.pdf) for the detail.
