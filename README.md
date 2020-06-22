## GISTIC2 for Canine CanFam3.1 genome

#### Changes to original code by

* [Emmanuel Martinez](https://github.com/jemartinezledes)
* [Samir B. Amin](https://github.com/sbamin)

For canine (canFam3.1), gistic2 run needs a few changes. You will need matlab to change code and recompile top level executable module or script, `gp_gistic2_from_seg`. I am using **v2.0.22** from ftp://ftp.broadinstitute.org/pub/GISTIC2.0/

**PS:** If you are unable to re-compile (so as we) matlab script using updated GISTIC2 version, try using pre-compiled matlab script, `gp_gistic2_from_seg_upd`. For this to work, you will need to use *canfam3_1_order.mat* and may require to downgrade Matlab runtime to MCR/8.0 and gistic2 to [v2.0.22](ftp://ftp.broadinstitute.org/pub/GISTIC2.0/all_versions/). **Please note** that our pre-compiled matlab script is based on relatively older reference genome annotations and cytoband information for CanFam 3.1 (~2011). Please make sure to **manually verify GISTIC2 peaks** by loading input copy-number segments into genome browser like [IGV](http://software.broadinstitute.org/software/igv/home) or [JBrowse](https://jbrowse.org/).

### Recompile `gp_gistic2_from_seg`

Following are a few changes I made after getting input from a colleague, [Emmanuel Martinez](https://github.com/jemartinezledes).

*   Make a **refgene.file** in `.mat` format containing three objects You can check valid format by  opening *hg19.mat* in matlab. 

1. cyto (1x40): containing cytoband. Since canFam has basic cytoband, it will not have detailed (1p11.1, 1p22, etc.) band info as with hg19.mat. canine cytoband can be fetched from UCSC table browser or ftp. https://link.sbamin.com/2ztpbWM or http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/cytoBandIdeo.txt.gz

2. rg (1x<number_of_genes>): gene metadata based on canFam refGene/Ensemble gtf file

3. rg_info (1x1): is additional metadata

Save these three objects as `canfam3.1.mat`

*   Now, extract GITIC2 v2.0.22 tarball. Structure would be something similar as noted here: ftp://ftp.broadinstitute.org/pub/GISTIC2.0/README.txt Now, we edit `source/RefGeneInfo.m` file, line 15 to change `nchr` from 24 to 40, and at line 29, add corresponding canine chromosomes c(25:38, 'MT', 'X').

```
RGI.chr.symb = {'1','2','3','4','5','6','7','8','9','10',...
                '11','12','13','14','15','16','17','18',...
                '19','20','21','22','23','24','25','26',...
                '27','28','29','30','31','32','33','34',...
                '35','36','37','38','MT','X'};
```

>Please make sure to have chr sequence order identical to one in cyto object, which was specified under `canfam3.1.mat`.  

*   Also, edit line 42 to 44 to match number of chromosomes: This may vary whether or not you like to add/remove autosomes, sex and MT chr.

```
RGI.txt2num('39') = 39;
RGI.txt2num('40') = 40;
RGI.chr.autosomal = (1:40)<39;
```

*   That's it! Recompile `gp_gistic2_from_seg` module which is the main executable under gistic2 shell wrapper.
  * [Here is guide](https://www.mathworks.com/help/compiler/create-and-install-a-standalone-application-from-matlab-code.html) to package matlab application. For the first step under an option, *add main file*, provide path to `gp_gistic2_from_seg.m`, add relevant details, keep LICENSE unchanged, i.e., it is copyrighted software by Mermel et al., Broad Institute, Cambridge, MA, USA. 
  *  I will push gistic2 code for canine sometime soon at https://github.com/sbamin Just fixing some of canfam3.1.mat issues to better annotate gistic plots. Feel free to comment/edit here for improvements/bugs.

### Segment file

*   View [make_segments.R](make_segments.R) for example script.

*   I use following script to generate segment file. It is well-defined at ftp://ftp.broadinstitute.org/pub/GISTIC2.0/

```r
## read each segment file, convert non-integer chromosome names to integers.

## create marker column based on the original window size used for making mappability bigwig file, i.e.,
## generateMap.pl -w 100 -i bowtie_index/CanFam3_1.fa "${REF_FASTA}" -o bigwigs/CanFam3_1.map.ws100bp.bw
```

~~convert log2 cn to log2-1 cn for GISTIC, so as log2(2) -1 == 0~~

**Note:** HMMcopy derived values are already **log2 copy ratio**, i.e., `log2(T/N) = log2(T) - log2(N) = log2(T) - 1`. So, it does not require substraction by 1, as I wrongly stated earlier. Read more on format at ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm and https://www.biostars.org/p/174382/#175590

```r
gistic2_seg_file = dplyr::bind_rows(lapply(1:length(segfiles), function(i) {
    aux_data <- read_tsv(segfiles[i], col_types = "cddid") %>%
                filter(!chr == "MT") %>%
                mutate(chri = match(chr, info_chr$X1),
                       start = as.integer(start),
                       end = as.integer(end),
                       markers = as.integer((end - start + 1)/100),
                       segvalue = median,
                       Sample = sampleids[i]) %>%
                dplyr::arrange(chri) %>%
                dplyr::select(one_of(c("Sample", "chr", "start", "end", 
                                "markers", "segvalue")))
    return(aux_data)
    }))

head(gistic2_seg_file)
colnames(gistic2_seg_file) <- c("Sample", "Chromosome", "Start Position", "End Position", "Num markers", "Seg.CN")

write_tsv(gistic2_seg_file, path = sprintf("%s/%s.tsv", getwd(), outfile))
```

### Marker file

*   marker file is not required since GISTIC v2.023 but we are using v2.022
*   Following [a post by Ming Tang](http://crazyhottommy.blogspot.com/2017/11/run-gistic2-with-sequenza-segmentation.html) to get marker file.

```sh
## marker file:

sed '1d' gistic2_segments.tsv | cut -f2,3 > markers.txt
sed '1d' gistic2_segments.tsv | cut -f2,4 >> markers.txt

## sort the files by chromosome, take the unique ones and number the markers.

cat markers.txt | sort -V -k1,1 -k2,2nr | uniq | nl > markers_gistic.txt
rm markers.txt
```

### Example run

View [gistic_run.sh](gistic_run.sh) file for details.

Samir
