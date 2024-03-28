# Accelerating K-Means Clustering for Documents<br> with an Architecture-Friendly Pruning Method
This repository contains supplemental materials including an [additional document](./supp.pdf) 
and a code set for appying our K-means clustering algorithm, AF-ICP, 
to **Large-scale and High-dimensional sparse data sets** such as 
the 8.2M-sized PubMed data set and comparing it with the other algorithms, 
ICP, TA-ICP, and CS-ICP.
The codes are implemented with C.

## Requirements for executing codes
1. OS: CentOS 7.6 and later
2. g++ (GCC): >= 8.2.0
3. perl: >= 5.16
4. perf: 3.10
5. bzip2 (optional)

## Quick start: AF-ICP in five iterations
1. Prepare the 8.2M-sized PubMed data set with a procedure in [dataset](./dataset).<br>
This procedure creates ``./dataset/pubmed.8_2M.db`` 
that is avilable for the codes in this repository.
You can download [pubmed.8_2M.db.bz2](http://prec4306.kanagawa-u.ac.jp/hp/pubmed.8_2M.db.bz2)
if you fail to download the original data (docword.pubmed.txt) from UCI machine learning repository.
Then, execute ``bzip2 -d pubmed.8_2M.db.bz2`` to extract the ``pubmed.8_2M.db``.
2. Execute ``make -f Makefile_itr5_aficp`` in ``./src``.<br>
 This makes ``./bin/itr5_aficp`` object in your system.
3. Execute the perl script ``./itr5_exeAFICP_8.2Mpubmed_perf.pl`` in ``./exe``.<br>
 The 8.2M-sized PubMed data set is loaded from ``./dataset/pubmed.8_2M.db`` (3.8GB)
 in around two minutes
 and given K=10,000, AF-ICP is executed with 50-thread parallel processing.<br>
 A log file is generated in ``./Log``.<br>

## Compare AF-ICP with other algorithms, ICP, TA-ICP, and CS-ICP
Go to [Comparison](./Comparison).

## License
Please check [LICENSE](LICENSE_v2.1.pdf) for the detail.
