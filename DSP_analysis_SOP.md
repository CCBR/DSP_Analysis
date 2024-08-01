# DSP Analysis QC SOP

---


## Index of actions for the QC Report

### 1) Check Inputs and Setup

Format input files and creating a the GeoMx object. This step will map probe IDs to gene names and combine count and metadata files

### 2) Run QC and Filtering

Add QC flags for specific QC parameters and filter AOIs and probes based on those flags. There will be a table for both AOIs and probes to review. 

### 3) Run Normalization and Check for Technical Error

Apply a normalization method to the raw count data for removal of technical error. Check the results of normalization using PCA, MA, or Concordance plots

### 4) Generate the QC Report

Use the Quarto sheet to generate a HTML report of all initial QC findings including Sankey Plot, QC, Filtering, and Normalization

## Index of actions for the Post-QC Analysis

### 1) Run Differential Expression

Differential expression analysis to compare gene counts between groups of AOIs and idenity differentially expression genes

### 2) Create Plots 

A number of different plots can be created based on the annotations provided. From differential expression a volcano plot can be generated for each contrast. Boxplots and heatmaps can also be generated for all AOIs or subsets based on annotations. 

### 3) Run Pathway Analysis

Gene set enrichement analysis for identifying pathways enriched in the differnetiall expressed genes of a particular differential expression comparison

---
---

#### Key Terms

**ROI**
: Region of Interest. 

>This is the selected area on a slide that AOIs are captured from. Each ROI will be assigned a number that us unique within a slide but NOT unique if there are multiple slides (each slide starts labeling ROIs from "001")

**AOI**
: Area of Interest. 

>The area within an ROI where probes are detached/collected from using UV light. Each line of the annotation file represents one AOI. Sometimes referred to as a "segment"


---
---

## 1. Inputs and Setup

### Action Step: Create the GeoMxSet Object

Run the function `studyDesign` from the DSPWorkflow package to create a GeoMxSet object. The GeoMxSet object will combine the metadata and read counts for easier down-stream processing.

Example:

`
sdesign.list <- studyDesign(dcc.files = dcc.files, 
                                pkc.files = pkc.files,
                                pheno.data.file = pheno.data.file,
                                pheno.data.sheet = "annotation",
                                pheno.data.dcc.col.name = "Sample_ID",
                                protocol.data.col.names = c("ROI"),
                                experiment.data.col.names = c("panel"),
                                slide.name.col = "slide name", 
                                class.col = "class", 
                                region.col = "Region", 
                                segment.col = "segment",
                                area.col = "area",
                                nuclei.col = "nuclei", 
                                sankey.exclude.slide = FALSE, 
                                segment.id.length = 10)
`

---

### Input Files Information and Requirements

---

#### The Annotation File

---

An excel file (.xlsx) where each row provides metadata for a specific AOI. The file needs the following required fields:

**class**

>The largest grouping for samples.
>
>Example: “Cancer” or “Normal”

**region**

>The anatomical region.
>
>Example: “TME” or “Tumor” or “Vessel”

**segment**

>The mask applied
>
>Example: “Full ROI” or “PANCK”

**slide_name**

>The name of the slide used in the DSP experiment.

**SampleID**

>The name of the sample which corresponds to the .dcc file name for that sample.
>
>Example: "DSP-1001660021204-D-A02"

**panel**

>The probe set used. It does not need to be the same as the PKC file name.
>
>Example: (v1.0) Mouse NGS Protein Core

#### Additional (optional) fields:

**nuclei**

>The nuclei count per AOI

**area**

>The area count per AOI

**Additional annotations**

>These are fields which contain data that can be used to create groupings for things like differential expression analysis or PCA plots


*Note on the annotation file*

The row order matters in this file when considering the No Template Control (NTC) wells. Each plate will typically contain 1 NTC well. Each NTC well should be listed before the rest of the AOIs for a given plate. Example:

| AOI | Plate |
| ----------- | ----------- |
| No template control | plate 1 |
| AOI 1.1 | plate 1 |
| AOI 1.2 | plate 1 |
| No template control | plate 2 |
| AOI 2.1 | plate 2 |
| AOI 2.2 | plate 2 |


*You will not have these columns in your annotation file, it simply shows the row order and how to arrage the NTCs*

---

#### DCC files

---

DCC files (extension ".dcc") are count files that contain probe IDs, counts for each probe ID, and sequencing QC information.

There will be one file per AOI captured, as well as files for the No Template Controls (NTCs).

Here is an example of a .dcc file:

```<Header>
FileVersion,0.02
SoftwareVersion,"GeoMx_NGS_Pipeline_2.3.3.10"
Date,2023-11-27
</Header>

<Scan_Attributes>
ID,DSP-1001660021203-C-A02
Plate_ID,1001660021203
Well,A02
</Scan_Attributes>

<NGS_Processing_Attributes>
SeqSetId,LH00236:67:22FLTKLT3
Raw,93081060
Trimmed,92645410
Stitched,90944511
Aligned,88961700
umiQ30,0.9923
rtsQ30,0.9961
</NGS_Processing_Attributes>

<Code_Summary>
RTS0020777,9288
RTS0020785,14017
RTS0020787,1141367
RTS0020790,43326
RTS0020796,2197
```

*Note on .dcc files*

Above you can see the .dcc provides QC information as well as counts for each probe, denoted in the "Code_Summary" section as an RTS number.

---

#### PKC files

---

A mapping file for probe ID to gene name. There will be one PKC file for each probe set used. If multiple probe sets were used there will be multiple PKC files.

PKC files are available for download from the [Nanostring website](https://nanostring.com/products/geomx-digital-spatial-profiler/geomx-dsp-configuration-files/)

---

## 2. QC and Filtering

### Action Step: Run QC Preprocessing

Run the `qcProc` function from the DSP Workflow package.

Example:

`
qc.output <-  qcProc(object = sdesign.list$object,
                        min.segment.reads = 1000, 
                        percent.trimmed = 80,    
                        percent.stitched = 80,   
                        percent.aligned = 80,    
                        percent.saturation = 50, 
                        min.negative.count = 3,   
                        max.ntc.count = 1000,     
                        min.nuclei = 200,         
                        min.area = 1000,
                        print.plots = TRUE)
`

The output will be a GeoMxSet object as well as tables for flagged segments and probes. 

*Note: AOIs will be filtered out based on the cutoffs defined in the QC parameters*

---

### QC Parameters' Definitions

---
 
#### Minimum Segment Reads

>The number of raw read per segment (default - 1000 reads)

#### Minimum Trimmed Reads

> The percentage of all reads in a segment that remain after removing adapter sequences (default - 80 percent)

#### Minimum Stitched Reads

> The percentage of all reads in a segment that remain after combining paired end reads that have overlapping sequences (default - 80 percent)

#### Minimum Aligned Reads

> The percentage of all reads in a segment that remain after aligning the stitched reads with the RTS ID barcode from the reference assay (default - 80 percent). 
>
>After alignment, deduplication occurs by removing PCR duplicates to create deduplicated reads 

#### Minimum Saturation

>Sets the minimum percent of sequencing saturation allowed. The percent of sequencing saturation is calculated as: 
>
>1 - (deduplicated reads/aligned reads) x 100
>
>100% sequencing saturation indicates a representative sample, while 0% sequencing saturation indicates that all reads were
unique. Values below 50% may need to be resequenced (default - 50 %)

#### Minimun Negative Count

>The minimum number of reads per segment from the geometric mean of negative probes (default - 10)

#### Max NTC Count

>The maximum number of reads in the No Template Control (NTC) well of the 96 well plate that was sequenced. NTC wells contain no sample and thus will indicate if there is contamination in the library prep (default - 1000)

#### Minimum Nuclei

>The minimum number of nuclei counted for a segment (default 200)

#### Minimum Area

>The minimum area of a segment (default 16000)

---

### Action Step: Run Filtering

Run the Filtering section of the Quarto document. This step is not currently available as a single function.

Filtering will first occur for AOIs (also called segments) followed  by filtering for probes (which respresent genes or negative controls)

The output will be a GeoMxSet object with selected AOIs and probes removed (determined by the QC and Filtering step). There will also be a series of graphs and tables on the detection rate of probes for each AOI and for specific annotation groups.

## 3. Normalization

### Action Step: Run Normalization

Run the `geomxNorm` function from the DSPWorkflow package. You will select a type of normalization to apply, either quartile-3 (q3) or negative (neg). The best practice is to run both normalization types for comparison.

Example:

```
 q3.normalization.output <- geomxNorm(
                                  object = object.gene.filtered, 
                                  norm = "q3")
    
    
neg.normalization.output <- geomxNorm(
                                  object = object.gene.filtered, 
                                  norm = "neg")
```

The output will be an object with normalized read counts.


#### Normalization Information

https://rdrr.io/github/Nanostring-Biostats/GeomxTools/src/R/NanoStringGeoMxSet-normalize.R

## 4. Post-Normalization Visualization

### Action Step: Compare Normalized Reads to Background

The `geomxNorm` function will output several plots that compares the Q3 value of each AOI to the mean counts of the negative probes, referred to as the negative background. Evaluate these plots for AOIs that have Q3 values that do not rise above the negative background, or do not follow the same trend of other AOIs of distance from Q3 to background. These plots are labeled by the annottion field region to check for bias. 

### Action Step: Create PCA Plots for Normalization Types

Run the `pca` function to generate principle component (PC) data for each normalization type. Run `biplot` to plot the two highest PCs as a dot plot where each dot represents an AOI and is labeled using one annotation field. Following the QC report, the four main annotaiton fields of slide_name, class, region, and segment are used to create one plot per normalization type per annotaiton field for 8 total PCA plots. Also create PCA plots any additional annotaiton fields of interest. 

Compare the normalization types for each annotation field. Depending on the experiment, the AOIs should cluster based on the annotation fields that represent biological variation and not technical variation. 

For example, AOIs clustering based on the anatomical region can be explained by gene expression differences. AOIs clustering based on the slide may represent differences in the slide preparation that affect the probe counts, which is what normalization aims to remove before down stream analysis.


### Action Step: Create MA Plots

Run the `make_MA` function (sourced from the file `DSP_QC_functions.R`) using two annotation groups to compare. Usually this comparison will be meaningful to the experimental design, such as region A versus region B. Compare the Q3 and Negative normalization MA plots to the raw counts to evalue the effect of normalization on the distribution. 

Example:

```
MA.plots.q3 <- make_MA(contrast.field = contrast.field, 
                       condition.label = condition.label, 
                       reference.label = reference.label, 
                       log.counts = log.counts, 
                       raw.log.counts = raw.log.counts, 
                       annotation = annotation.MA)

grid.draw(MA.plots.q3)
```

In the Quarto document, you will need to define the annotations to use in the construction of the MA plots. See the Quarto report example for more information.


