# Experimenal Notes

This document is to keep track of decisions as and when they are made. Format should be data, question, decision, and explanation/justification for why the decision was made

## 05-09-2022

#### How should we filter the data?

We first want to filter outcomes that are noise. We want to keep outcomes that:
- Feature commonly across a number of gRNAs. 
- Constitute a significant part of the repair outcome profiles of any gRNA
- We might want to filter out genes which have <3 gRNAs

- We would like some measure/proxy for the quality of a particular gRNA sequence. Right now, we are using the number of reads.

- We need to filter barcodes such that the barcodes in use are common to all sample sequences

- Our aim should be to filter out RARE outcomes. Not keep "interesting" outcomes. 

We need to examine each sample seperately. 




