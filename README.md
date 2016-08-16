# Project : Compare LEON and GZIP for the FastQ compression

__Table of Contents__

- [Project : Compare LEON and GZIP for the FastQ compression](#)
	- [Just a little story](#just-a-little-story)
	- [What is currently done?](#what-is-currently-done)
		- [GZIP](#gzip)
	- [What is LEON ?](#what-is-leon-)
	- [Comparison between compression of fastq by Gzip and LEON](#comparison-between-compression-of-fastq-by-gzip-and-leon)
		- [Compression ratio fo each tools](#compression-ratio-fo-each-tools)
		- [Compression ratio depending on the size of FASTQ](#compression-ratio-depending-on-the-size-of-fastq)
		- [Time of compression/decompression](#time-of-compression-decompression)
	- [Is LEON compression have an impact on SNPs/Indels calling ?](#is-leon-compression-have-an-impact-on-snpsindels-calling-)
	- [Related publications](#related-publications)

## Just a little story

In the 1970's Sanger, Maxam, Gilbert and colleagues developed a rapid method to sequence the DNA. 
Twenty years after, sequencing by Sanger method is the more common way and allows the first whole genome sequencing for _Haemophilus influenzae_, in 1995.
In 2004, almost thirty years after Sanger has developed his method, the Human Genome Project sequenced the first whole Human genome. Since this step, sequencing methods have changed and the Next Generation Sequencing (NGS) have emerged.
In 10 years the cost and the time to sequence a whole genome considerably decrease. NGS technologies allow to sequence routinely a large number of samples.
So, the amount of data generated by NGS substantially increase during the last decade and the storage and transmission of these data are a major concern for now.

![Graph from SRA (http://www.ncbi.nlm.nih.gov/Traces/sra/) 2016-08-08](https://github.com/Char-Al/bench_leon/blob/master/images/NGS_data.png "The SRA database, wich contains a large part of the world wide sequencing, is growing very fast and now contains almost 6 petabases (date : 2016-08-08)")

## What is currently done?

### GZIP

Now, the common way to compress those data is the GZIP format. GZIP is based on the Deflate algorithm, actually it is the combination of the Huffman coding and the LZ77 algorithm ([more explanation here](http://www.zlib.net/feldspar.html)).
This algorithm have been developed in order to compress text data, which means data with a large set of characters.

## What is LEON ?

LEON is a new software to compress data issue from NGS (Fasta and FastQ).
LEON shares similarities with approaches using a reference genome to compress files.
The LEON approach build the reference _de novo_, contrary to the other algorithms, with a _de Bruijn Graph_ where the pieces are _k-mers_.
The _de Bruijn Graph_ is heavy and it has to be stored with the compressed data, so the size could be a problem.
To figure out this problem, the _de Bruijn Graph_ needs a good parametrization and the implementation is based on probabilistic data structure in order to reduce its size. Based on Bloom filters the _de Bruijn Graph_ is efficient to store large data.

![LEON method overview (from : Reference-free compression of high throughput sequencing data with a probabilistic de Bruijn graph)](https://github.com/Char-Al/bench_leon/blob/master/images/LEON_overview.png "LEON method overview (from : Reference-free compression of high throughput sequencing data with a probabilistic de Bruijn graph)")

## Comparison between compression of fastq by Gzip and LEON

With this little magic script, we produce some awesome graphs to compare the efficiency of GZIP and LEON.
To compare these two softwares, we are interested in the global rate of compression, the rate of compression depending the size of the initial FastQ and the time of compression/decompression.
We use FastQ from Human data with size between 100 Mo and 26 Go.

### Compression ratio fo each tools

We can see that the ratio compression of LEON is better than GZIP, regardless the size of the FastQ.
In addition the _"lossy"_ mode of LEON have a ratio between 90 and 95% in each cases, almost 15% more than the other tools.
There is no significant differences between the level 6 and 9 of GZIP, but these two have a wider variations.

![Boxplot comparant les taux de compression de gzip et LEON avec différentes options](https://github.com/Char-Al/bench_leon/blob/master/example/boxplot_compression.png "Boxplot comparant les taux de compression de gzip et LEON avec différentes options")

### Compression ratio depending on the size of FASTQ

Now, we focus on the compression rate depending of the size of the original FastQ.
We can notice a peak at 18 Go, corresponding to FastQ files that have larger reads (125 vs 100 pb).
However this peak cannot change the analysis because all softwares show the same effect.
Regarding these results, we can say that LEON is more efficient than GZIP, especially the _"lossy"_ mode which is very stable in each case.

![Evolution du taux de compression en fonction de la taille des fastQ d'origine](https://github.com/Char-Al/bench_leon/blob/master/example/point_compression.png "Evolution du taux de compression en fonction de la taille des fastQ d'origine")

### Time of compression/decompression

The time of compression and decompression depends of the size of initial file.
The LEON _"lossy"_ and _"lossless"_ mode and GZIP level 6 have similar time for the compression, while the time of GZIP level 9 for compression is longer (in some cases more than 2 times).
LEON is less efficient for decompression and all GZIP level have almost the same time for decompression.

![Evolution du temps de compression en fonction de la taille des fastQ d'origine](https://github.com/Char-Al/bench_leon/blob/master/example/point_time.png "Evolution du temps de compression en fonction de la taille des fastQ d'origine")

## Is LEON compression have an impact on SNPs/Indels calling ?

To studying the impact of compression we will used a set of 12 FastQ issue from Human data.
All FastQ have been compressed by LEON _"lossy"_, LEON _"lossless"_, and GZIP (default : level 6), then decompressed and recompressed with GZIP (default).
The last recompression is usefull to perform the variant calling with Nenufaar v2.3 (pipeline by David BAUX).

The next table show the number of SNV and indels call for each VCF issue from the three compression methods.
We notice that LEON _"lossy"_ mode have 23 SNV and indels which differ than others.
Indeed, 12 ponctual mutations are found only with GZIP and LEON _"lossless"_ mode and 11 with LEON _"lossy"_.
However most of these 23 mutations are in repeated sequence and this difference may be caused by a shift of few nucleotides.

VCF                   | PASS | OTHER |  NA  | TOTAL (PASS+OTHER)
--------------------- | :--: | :---: | :--: | :----------------:
__GZIP__              |  759 |  110  |  12  |         869
__LEON _"lossless"___ |  759 |  110  |  12  |         869
__LEON _"lossy"___    |  767 |  103  |  11  |         870

The following chart show the differences between the AB (Allelic Ballance) of each SNV/indel.
The AB from the VCF issue of GZIP (VCF1).
We notice than there is no differences between GZIP and LEON _"lossless"_ (VCF2).
With LEON "_lossy"_ mode there is some differences.
The AB of 84 indels, on 192 indels in totals, is different but most of these (66) have an AB wich differ less than 2%.
We can have the same conclusion with SNV with 344 wich differ, on a total of 666 SNV, and 318 wich differ less than 2%.

![Differences of variant calling on GZIP, LEON lossless and LEON lossy file](https://github.com/Char-Al/bench_leon/blob/master/example/callingDiff.png "Differences of variant calling on GZIP, LEON lossless and LEON lossy file")


## Related publications
* International Human Genome Sequencing Consortium. [Finishing the euchromatic sequence of the human genome](http://www.nature.com/nature/journal/v431/n7011/full/nature03001.html) Nature 431, 931–945. issn: 1476-4687 (Oct. 2004).
* Fleischmann, R. D. et al. [Whole-genome random sequencing and assembly of Haemophilus influenzae Rd](http://science.sciencemag.org/content/269/5223/496.long) Science (New York, N.Y.) 269, 496– 512. issn: 0036-8075 (July 1995).
* Sanger, F., Nicklen, S. & Coulson, A. R. [DNA sequencing with chain-terminating inhibitors](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC431765/) Proceedings of the National Academy of Sciences of the United States of America 74, 5463–5467. issn: 0027-8424 (Dec. 1977).
* Benoit, G. et al. [Reference-free compression of high throughput sequencing data with a probabilistic de Bruijn graph](http://www.biomedcentral.com/1471-2105/16/288) BMC bioinformatics 16, 288. issn: 1471-2105 (2015).
* Van Dijk, E. L., Auger, H., Jaszczyszyn, Y. & Thermes, C. [Ten years of next-generation sequencing technology](http://www.sciencedirect.com/science/article/pii/S0168952514001127) Trends in genetics: TIG 30, 418–426. issn: 0168-9525 (Sept. 2014).
