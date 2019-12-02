# Exploring 16S data from Oral treatments with probiotics and live Salmonella

* 11/30/2019
* local:/Users/severin/Desktop/Projects/Mellata/01_qiime

Now that we have removed lane effect, how do the heatmaps look and can we identify groups of OTUs that may explain some of the differences between treatment groups.


## Gneiss analysis of frequency table

We can explore how different fixed factors contribute to the abundances in the frequencies.a

### Gneiss Correlation-clustering

This clustering method uses the Ward hierarchical cluster.  If two microbes are highly correlated across all of the samples, then the distance will shrink close to zero

```
qiime gneiss correlation-clustering \
  --i-table table-dada2.qza \
  --o-clustering hierarchy.qza
```

### Building linear models using balances¶

The isometric log ratio (ILR) transform computes the log ratios between groups at each node in the tree.  These balances avoid situations where one very abudundant species can cause a comparison between organisms to appear different.

```
qiime gneiss ilr-hierarchical \
  --i-table table-dada2.qza \
  --i-tree hierarchy.qza \
  --o-balances balances.qza
```


## Building linear models using balances¶

I changed Probiotics to AProbiotics to make it alpha numerically the first of the treatment categories so that it uses it to compare the other treatments to.  OLS = Ordinary least squares

I also added a column to the sample-metadata.tsv file that clusters the samples by the unknown effect observed by the unifrac-emperor plot.

```
qiime gneiss ols-regression \
  --p-formula "Description1+Unknown" \
  --i-table balances.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization regression_summary_lane.qzv

  qiime tools view regression_summary_lane.qzv
```

### Output

We can see that the model fits this data very well. The live salmonella treatments to probiotic account for ~20% of the variation and the unknown factor accounts for about 40% of the variation.  In addition the predicted values (red) cluster with the raw data (blue) and the residuals are small (<6) compared to the variation (~25)

![SimpleLinearRegressionSummary](/Notebook_Severin/02_gneissAnalysis/SimpleLinearRegressionSummary.png)

![projected prediction](/Notebook_Severin/02_gneissAnalysis/Projectedprediction.png)

There may be some concern about the large unknown variation but as we saw in the emperor plot, the treatment affect is orthoganol to the unknown effect.  This is further shown in the significant fdr corrected p-value balances.  

There are only 5 balances that are fdr corrected p-value significant for the unknown factor. Of these only y290 just barely passes the 0.05 threshold for sinificance for T.Probiotics-AND-LiveSalmonella compared to probiotics.  I state this so we know that we can be confident in the groups of OTUs we explore later on that distinguish the treatement groups. We will know that it is real and not influenced by the unknown factor.


```
      Intercept	Description1[T.LiveSalmonella]	Description1[T.Probiotics-AND-LiveSalmonella]	Unknown[T.B]
y257	1.27E-06	0.705070766	0.128992934	1.42E-07
y290	3.51E-13	0.342813711	0.049940587	0.01665184
y459	1.59E-06	0.741851534	0.53397014	0.021635915
y81	0.295894629	0.866739752	0.883206209	0.029062446
y15	0.090943511	0.449361564	0.54632248	0.029370752
y21	0.287344409	0.038133526	0.273428412	0.029742541
```



## Heatmap Treatments

There are several balances that appear to differentiate the treatment groups based on intensity on the heatmap and we can confirm based on the fdr-corrected pvalues.

```
qiime gneiss dendrogram-heatmap \
  --i-table table-dada2.qza \
  --i-tree hierarchy.qza \
  --p-ndim 30 \
  --m-metadata-file sample-metadata.tsv  \
  --m-metadata-column Description1 \
  --p-color-map seismic \
  --o-visualization heatmapTreatment.qzv
  qiime tools view heatmapTreatment.qzv
```
![](/Notebook_Severin/02_gneissAnalysis/HeatmapTreatment.png)
## Heatmap Unknown

unlike the treatement heatmap, the unknown heatmap appears to be random.  This coinsides with the very few significant balances.

```
qiime gneiss dendrogram-heatmap \
  --i-table table-dada2.qza \
  --i-tree hierarchy.qza \
  --p-ndim 30 \
  --m-metadata-file sample-metadata.tsv  \
  --m-metadata-column Unknown\
  --p-color-map seismic \
  --o-visualization heatmapUnknown.qzv
  qiime tools view heatmapUnknown.qzv
  ```

![](/Notebook_Severin/02_gneissAnalysis/HeatmapUnknown.png)


## Significant balances and OTUs

Let's focus on only those balances with an fdr corrected pvalue of less than 1e-5

### Description1[T.LiveSalmonella] compared to Probiotic
```
Intercept	Description1[T.LiveSalmonella]	Description1[T.Probiotics-AND-LiveSalmonella]	Unknown[T.B]
y6	5.64E-05	2.09E-15	2.09E-15	0.625753313
y0	8.42E-15	1.65E-13	2.75E-14	0.591743441
y27	7.44E-10	2.56E-08	6.76E-08	0.460321344
y2	0.285657342	8.11E-08	0.173854082	0.449523148
y14	0.533871326	5.83E-07	1.36E-15	0.382380621
y207	1.08E-05	2.03E-05	2.03E-05	0.30990384
y28	0.076478361	2.55E-05	0.060795012	0.566888476
y394	4.54E-10	7.26E-05	0.057777827	0.44022758
y92	0.280535721	9.97E-05	0.191208525	0.137629888
```

Only y2, y28,y394, and y92 separate the two treatment groups based on a 0.05 corrected pvalue significance. With y2 and y92 being very significant in  T.LiveSalmonella and not even close to significance in T.Probiotics-AND-LiveSalmonella.

### OTUs that distinguish T.LiveSalmonella

#### y2


#### y92





### Description1[T.Probiotics-AND-LiveSalmonella] compared to Probiotic

```
y14	0.533871326	5.83E-07	1.36E-15	0.382380621
y6	5.64E-05	2.09E-15	2.09E-15	0.625753313
y0	8.42E-15	1.65E-13	2.75E-14	0.591743441
y27	7.44E-10	2.56E-08	6.76E-08	0.460321344
y207	1.08E-05	2.03E-05	2.03E-05	0.30990384
y174	3.29E-05	0.000102031	3.91E-05	0.746212239
y294	3.94E-18	0.000887142	9.48E-05	0.448059245
```

Nothing really distinguishes T.Probiotics-AND-LiveSalmonella as having a unique set of OTUs as everything in at 1e-5, however if we go down a little farther we do find a the next two balances where it is significant in T.Probiotics-AND-LiveSalmonella and not significant in T.LiveSalmonella.

```
y43	0.000260666	0.266408348	0.000320472	0.145755106
y447	1.65E-17	0.836744138	0.00099771	0.836744138
```

### OTUs that distinguish T.Probiotics-AND-LiveSalmonella

#### y43


#### y447



### OTUs that distinguish probiotic from livesalmonella and livesalmonella+probiotics

#### y0

#### y6

#### y14

#### y27
