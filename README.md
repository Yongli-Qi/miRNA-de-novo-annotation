# miRNA-de-novo-annotation

SYNOPSIS: 
    Establishment of miRNA library after using ShortStack program (Johnson NR, Yeoh JM, Coruh C, Axtell MJ. (2016). G3 6:2103-2111.
    doi:10.1534/g3.116.030452)

USAGE: 
    python miRNA_de_novo_annotation.py [three-letter abbreviation of species name] [reads]

Core OUTPUT FILE:
    info6.txt 
    
      The file info6.txt is a plain-text tab-delimited file that contains the core results of the analysis.
      The meaning of these columns is described below:
      
      1. Original loci: Locus in format Chr:Start-Stop, which is annotated from ShortStack program
      2. MiRNA name
      3. Chromosome
      4. Hairpin start	
      5. Hairpin end	
      6. Hairpin sequence	
      7. Hairpin structure	
      8. Strand	
      9. MiRNA start	
      10. MiRNA end	
      11. MiRNA sequence	
      12. MiRNA* start	
      13. MiRNA* end	
      14. MiRNA sequence	
      15. MiRNA* detected?	
      The NO in this column means the locus met the criteria of N15 in ShortStack (Passed all tests except 
      that the miRNA-star was not sequenced) and YES in this column met the criteria of Y in ShortStack (Passed 
      all tests including sequencing of the exact miRNA-star)
      16. TPTM: Total number of primary alignments normalized to reads per 10 million
      17. Location: 3'-end or 5'-end of hairpin
