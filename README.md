# Extracting microbial signal from host-dominated  metatranscriptomes

Antonin COLAJANNI <sup>1,2</sup>, Raluca URICARU<sup>2</sup>, Rodolphe THIÉBAUT<sup>1</sup>, and Patricia THEBAULT<sup>2</sup> 

<sup>1</sup> Univ. Bordeaux, INSERM, INRIA, BPH, U1219, F-33000 Bordeaux, France 

<sup>2</sup> Univ. Bordeaux, CNRS, Bordeaux INP, LaBRI, UMR 5800, F-33400 Talence, France

Corresponding Author: antonin.colajanni@u-bordeaux.fr 

**Keywords**

Metatranscriptomic, Metagenomics, Host-dominated metatranscriptomes, Microbial translocation


**Abstract** 

Human RNA-seq data originally generated for human transcriptome profiling are overwhelmingly dominated by host sequences, yet they often contain a small fraction of non-human reads that can be exploited for microbial detection. When such datasets are repurposed for secondary microbiome-oriented analyses, extracting and accurately classifying this weak microbial signal becomes technically challenging and no ready-to-use pipeline currently exists. 
In this study  we evaluate computational strategies for filtering host reads and classifying microbial transcripts in host-dominated RNA sequencing data. We compare assembly-based approaches similar to those used in a previous study focusing on microbial translocation, with state-of-the-art  assembly-free methods, and assess their respective strengths and limitations using simulated datasets reflecting low microbial abundance. Our results show that assembly-based methods yield accurate taxonomic predictions but struggle at low read depth, whereas assembly-free methods are  more robust in sparse settings at the cost of reduced precision.
To leverage the complementarity of both approaches, we propose a hybrid pipeline that integrates assembly-based and assembly-free classification. On simulated data, this hybrid strategy improves microbial classification performance compared to either approach alone. Application to a real human metatranscriptomic dataset analysed in a microbial translocation context illustrates the broader microbial signal captured by the hybrid approach, despite intrinsic challenges related  to the absence of reliable ground truth and to the risk of host read misclassification. 
Overall, our work provides a framework for extracting microbial signal from host-dominated human metatranscriptomes, enabling the reuse of existing transcriptomic datasets for microbiome-related analyses, including but not limited to microbial translocation studies.
