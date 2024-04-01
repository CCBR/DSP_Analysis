# DSP Analysis Documentation: QC Parameters

## DSP Analysis Video

https://university.nanostring.com/geomx-academy-ngs-data-analysis-for-rna/1029434


## QC Parameters Summary

### Nanostring Documentation

Nanostring University
- https://university.nanostring.com/page/document-library

DSP NGS Readout
- https://university.nanostring.com/geomx-dsp-ngs-readout-user-manual/1193408

DSP Data Analysis
- https://university.nanostring.com/geomx-dsp-data-analysis-user-manual/1163670




### QC Parameters
 
#### Minimum Segment Reads

The number of raw read per segment (default - 1000 reads)

#### Minimum Trimmed Reads

The percentage of all reads in a segment that remain after removing adapter sequences (default - 80 percent)

#### Minimum Stitched Reads

The percentage of all reads in a segment that remain after combining paired end reads that have overlapping sequences (default - 80 percent)

#### Minimum Aligned Reads

The percentage of all reads in a segment that remain after aligning the stitched reads with the RTS ID barcode from the reference assay (default - 80 percent). 

After alignment, deduplication occurs by removing PCR duplicates to create deduplicated reads 

#### Minimum Saturation

Sets the minimum percent of sequencing saturation allowed. The percent of sequencing saturation is calculated as: 

1 - (deduplicated reads/aligned reads) x 100

100% sequencing saturation indicates a representative sample, while 0% sequencing saturation indicates that all reads were
unique. Values below 50% may need to be resequenced (default - 50 %)

#### Minimun Negative Count

The minimum number of reads per segment from the geometric mean of negative probes (default - 10)

#### Max NTC Count

The maximum number of reads in the No Template Control (NTC) well of the 96 well plate that was sequenced. NTC wells contain no sample and thus will indicate if there is contamination in the library prep (default - 1000)

#### Minimum Nuclei

The minimum number of nuclei counted for a segment (default 200)

#### Minimum Area

The minimum area of a segment (default 16000)




