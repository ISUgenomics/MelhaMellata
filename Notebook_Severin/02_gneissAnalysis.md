# Exploring 16S data from Oral treatments with probiotics and live Salmonella vaccine

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
export YY='y2'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Description1 \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv

Numerator had 442 OTUs and is not likely the group that is distinguishing.

Denominator

```
Feature ID	0	1	2	3	4	5	6
0f96ba1e69e26e843b1da15f4393ce16	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
3132e062f465bfbc042219075df3caa1	k__Bacteria	 p__Firmicutes	 c__Erysipelotrichi	 o__Erysipelotrichales	 f__Erysipelotrichaceae	 g__Coprobacillus	 s__
1c5f74d5cbbe5b3aed73436377cc8cb2	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales
7b28c20e72c6c95b3e604f0849245770	k__Bacteria	 p__Verrucomicrobia	 c__Verrucomicrobiae	 o__Verrucomicrobiales	 f__Verrucomicrobiaceae	 g__Akkermansia	 s__muciniphila
5db6f4f8ef0dcb96ad17d5156b2f9d0e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
47b098c7edc1e61aef6d211df51b8d96	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__[Ruminococcus]	 s__
1fc9ec3a3fdc7d611ed88ee82f5c69f0	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Blautia	 s__
b7453d8e94f5556929a479a78cc538b9	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
bd5ab72117a5582f67a888c1f7e57ca7	k__Bacteria	 p__Firmicutes	 c__Erysipelotrichi	 o__Erysipelotrichales	 f__Erysipelotrichaceae	 f__Erysipelotrichaceae	 f__Erysipelotrichaceae
d87d635850bba0b65279b91092dec469	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
e935f4346c341bf7a59a3efb7dde55d5	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Dorea	 s__
664d7ec2732a1ad4ddad8a47483facec	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
0a6b9615944301a29e20464d2ebb5ae1	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
0c8c712659cf39aec314278ad552c0ed	k__Bacteria	 p__Firmicutes	 c__Erysipelotrichi	 o__Erysipelotrichales	 f__Erysipelotrichaceae	 g__	 s__
a42d3b0995ba1eb28bf897a341a7ac65	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
4c2396ab8eabb16c11d85f99ff4270f1	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
ea7823427b5854bc784da136aa27d777	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
26d9472c382472841670239be53bf945	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__[Ruminococcus]	 s__
45aea834beb1b1b53d8db337520b771e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
```

#### y92

```
export YY='y92'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Description1 \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv
```

This balance only had 1 OTU in each the numerator and the denominator so it isn't clear how useful it will be.

Numerator
```
0	1	2	3	4	5	6
e935f4346c341bf7a59a3efb7dde55d5	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Dorea	 s__
```

Denominator

```
0	1	2	3	4	5	6
47b098c7edc1e61aef6d211df51b8d96	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__[Ruminococcus]	 s__
```



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
```
export YY='y43'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Description1 \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv
```

Numerator
```
Feature ID	0	1	2	3	4	5	6
16ca5a38fd9fecf4e2581a98b6eb67f2	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
3120660c883750285d85c12f7d898f67	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
8f003929e210470eb85ae81509cb2dd2	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
bd5ed0bc219bedc17bec9cb7b6106f79	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
19442576619f36291c2c41b17009342e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Butyricicoccus	 s__pullicaecorum
cccbf8a72fe25e3311fd86cf66820a49	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
249a5e18a13e5f0567cdf7d9d32a577f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
cf42d6ef5b5f01ff0a0a2bf5b533682e	k__Bacteria	 p__Firmicutes	 c__Bacilli	 o__Bacillales	 f__Bacillaceae	 g__	 s__
c3f21101b8845fd1c974cb14578fb119	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Christensenellaceae	 g__	 s__
6923fabf0a9d16b03765a551acc70404	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
7ff0c8d6b1f2ab6ae6c70229e507bc18	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
47f770032218f91b352999ff1f663807	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Clostridium	 s__aldenense
77fa64d0d33043a8e44cd873c3cd4edb	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
33ddca6781100b91d670de9a543e5aa7	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
4412ce79670508864ba858560d7a5b42	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
8b414a059be1b8812e888b7c4d47bb2f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
5b84aff1766e777420629a2a961de7ea	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Butyricicoccus	 s__pullicaecorum
```

Denominator
```
Feature ID	0	1	2	3	4	5	6
2799ad5c90e307f0e897fecda5a2924e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
2cdac4ae837db323be3207da32af0649	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
a72bc4bbcb65e31aac628de5b9a9b0af	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Christensenellaceae	 g__	 s__
e5f3d35279d4d4f0c6836a68aef9d04f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
aa0d22e3cae06a85b97c6181d08c648a	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
5c0eb2cc0831a4b7304ec1fa20756159	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
78500a9f64f33ec04e40892a3f534176	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
6145d5800bec71cdd01246b13a463feb	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__[Ruminococcus]	 s__
ab1e4aeddde9eb4e80144c135051d941	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Eubacteriaceae	 g__Anaerofustis	 s__
0e226d197c2f253bcee5f15026512f90	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Dorea	 s__
c687a8d04aad1fff0bc7258d20d50db7	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
1a3af659e20f20ba69ce3ed6080723de	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
728c102c5aa9e0b548c64b2d6e0b1bf4	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
b83d7b512a42b855bf414420e91c6fd2	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
55e8fd9e2a3868ef69575fcb459fc7c5	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
9acb555d3c5dd77e9035a661e7ad1bb6	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Clostridiaceae	 g__Clostridium	 s__celatum
64c296ddb011fb2e7251ecadcb4a0148	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
a238176211bb5d089460b8d79136c127	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
1dc952db562f039bbc6a4dd86f4087e1	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
fbb2ccdf9eee6835b98ab9fad97f117a	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
50d3472ae91234c809fc3f45b8c26dfc	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Dorea	 s__
ea7064a0fdbb5228b7a0c9ed23e012d2	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
08c53fdf2c7b879ac5e38d5730b65fd9	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Christensenellaceae	 g__	 s__
cb3bfa7ac2b500c6a49d88ce9412dbea	k__Bacteria	 p__Firmicutes	 c__Erysipelotrichi	 o__Erysipelotrichales	 f__Erysipelotrichaceae	 g__cc_115	 s__
56eb9020ab3322ba96f782a4c25dddb9	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
```

#### y447

export YY='y447'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Description1 \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv

Numerator

```
Feature ID	0	1	2	3	4	5	6
2cd838736b8090ae5a6e5bc922b49999	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
e4ef1c7b616eb842c0ff7a327e83bd76	k__Bacteria	 p__Bacteroidetes	 c__Bacteroidia	 o__Bacteroidales	 f__Porphyromonadaceae	 g__Parabacteroides	 g__Parabacteroides
38cd1f968bf1b93f5e05ef1b5fb12ef4	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
233f1902fa0a8eb4c6a91321bc468a88	k__Bacteria	 p__Firmicutes	 c__Bacilli	 o__Lactobacillales	 f__Lactobacillaceae	 g__Lactobacillus	 s__hamsteri
23c8d53b5d6cb2f468b124365ab164a6	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
f13a69ad40cf466ff5bfb0c0609aafc6	k__Bacteria	 p__Actinobacteria	 c__Coriobacteriia	 o__Coriobacteriales	 f__Coriobacteriaceae	 g__	 s__
c0d5395792eadbf5f62e8ffb14fa0262	k__Bacteria	 p__Proteobacteria	 c__Alphaproteobacteria	 o__Rickettsiales	 f__mitochondria	 f__mitochondria	 f__mitochondria
99deb3c5ecb022ec05609ebd1112a557	k__Bacteria	 p__Bacteroidetes	 c__Bacteroidia	 o__Bacteroidales	 f__Bacteroidaceae	 g__Bacteroides	 s__
fde81acceeb6db492ac9f39724931f75	k__Bacteria	 p__Actinobacteria	 c__Coriobacteriia	 o__Coriobacteriales	 f__Coriobacteriaceae	 g__Eggerthella	 s__
2ef1e51ab1cf99a3c6417b05a060830e	k__Bacteria	 p__Firmicutes	 c__Bacilli	 o__Lactobacillales	 f__Lactobacillaceae	 g__Lactobacillus	 g__Lactobacillus
ba4925d004f3d04ff0ce5630dfe7e02d	k__Bacteria	 p__Bacteroidetes	 c__Bacteroidia	 o__Bacteroidales	 f__S24-7	 g__	 s__
1a4184a394e5f3db74169817b4d0b32f	k__Bacteria	 p__Firmicutes	 c__Bacilli	 o__Lactobacillales	 f__Lactobacillaceae	 g__Lactobacillus	 s__
c2b21f1d8e071fda68a09853a8cdcfdd	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Faecalibacterium	 s__prausnitzii
3bfd592733c45cdfe746a0840e332a3c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__eutactus
46308a5a9da9d9d102c8327b4d97ff62	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Faecalibacterium	 g__Faecalibacterium
61106a3c95b0b7dd9f989a76d12453b8	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
ad464b0f4f2cef292428ce4cca18f4dc	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
1f907a1c4b28365e0941c891d325ed29	k__Bacteria	 p__Bacteroidetes	 c__Bacteroidia	 o__Bacteroidales	 f__Prevotellaceae	 g__Prevotella	 s__
30735ab18567e115bd3f08fceaede146	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
5c60f84b3d868e8544e0b7805e20ae77	k__Bacteria	 p__Proteobacteria	 c__Gammaproteobacteria	 o__Pseudomonadales	 f__Moraxellaceae	 g__Acinetobacter	 s__
e5d61f7e72dcfaa67f72a58490bc6a0f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
c939dc55399371e18eb56a1b5ce84093	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales
ac414f95242bdb25d76ea1adbb11a1b3	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Blautia	 s__
84018a27b2b0ae81c662c50e204f6cf3	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
1bbb91863bac2d066d97bcb99dffd94d	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
```

Denominator

```
Feature ID	0	1	2	3	4	5	6
7c861f48eaca9c2cacdb0470b609b1d4	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
6b838da8ee8519ee41aee397dcce7d60	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae

```

### OTUs that distinguish probiotic from livesalmonella and livesalmonella+probiotics

#### y0 top two most significant balance to distinguish btween LS and probiotic treatment groups

```
export YY='y0'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Description1 \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv
```
Numerator has 460 OTUs

Denominator has 51 OTUs

```
Feature ID	0	1	2	3	4	5	6
77c2dc197e6b3dbebc4ee240c6a1c559	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
3eff8023f661e255741f87b115548289	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
6b5c5532185593c41966b2173ea9df1c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Clostridium	 s__lavalense
e657205bc51c1db7e6a38c33a084749d	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
5fbc483f79cd97d3a7253931f37c33e8	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
e9c788e08d5e698c0c90175a9770a1db	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Blautia	 s__
d46e2205f0c6ecf67b51f83d111c509c	k__Bacteria	 p__Proteobacteria	 c__Gammaproteobacteria	 o__Enterobacteriales	 f__Enterobacteriaceae	 f__Enterobacteriaceae	 f__Enterobacteriaceae
9cbf794b9c89435874fce22eba240120	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
816478dc0b3d7f8973f67307a0007320	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
1e8cd6c305d8f64ca2a2790405b64f6e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
9b6b08a24a54e5b482f5a45e03dcb30e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__[Ruminococcus]	 s__
3001460c772727615c37fd5f4c87489c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Faecalibacterium	 s__prausnitzii
100bb7f94035bd01a39239dd6803be6f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
ea6a2deed19ee7725916d60689c97b54	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
66b90e64bc0857f1f29cb63d373fcff0	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
2c3b179c770eca4dd9338a4d2c137202	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
708b000c48dbefb1454b79108d2a1f70	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
dc41692600de807e3e76f265cdc14075	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__[Ruminococcus]	 s__
99bccd58aaeb00f98e22add40be6e36b	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Peptostreptococcaceae	 f__Peptostreptococcaceae	 f__Peptostreptococcaceae
d6ceb660cd4c205d9389c35e067aac3f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
6b7a4e87dbafad40073befec06406e72	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Clostridiaceae	 g__Clostridium	 s__
8e9b12466534772cec6b7b991176c475	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
efb71540108125289bc7059eb3c1c43f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
bb667d3a920b293efb1dd114d721be64	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
7fc5ae42af9050001551d3dcd38bd57b	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Blautia	 s__
b2088770a003be218487d86b1d646231	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
8319d3eac517cdc6eac67429141a4179	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__[Ruminococcus]	 s__
4d3b09e575927f03bd03cfa6c30e1e58	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
f980a909308141ab17d7502f6c96d0ed	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
caa7ba533f97be6f54d95b2064cc33e5	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
49ff4b9ad904b7540926c72823626888	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
edef539dc6e5699b9f2c036009ad51ca	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
ac4274f46d9ec670be0cfb4126304a8d	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
970dd45e40cf8b847e2575fc1cb5a3eb	k__Bacteria	 p__Firmicutes	 c__Bacilli	 o__Lactobacillales	 f__Leuconostocaceae	 g__Weissella	 g__Weissella
eedf84efe5461c4221d334e78c8f6d82	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
c5619065fbb2332fe286b698e6b6c58d	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
c12636eadd51f39fc336f44d4cfac7f5	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Eubacteriaceae	 g__Anaerofustis	 s__
495fb0a927bfbbb88533b7127d146370	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
9c782dc640de8ecf83bef9fa1bf63ad0	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
23d4f160a1453f447e17669f62afbc86	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Blautia	 g__Blautia
0a39f269a92286a2fc485a177298da7a	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
b4f28264f150902936030ef29ce2a10a	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__[Mogibacteriaceae]	 g__	 s__
170a67d397000b76587d7756f3f6feb4	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
67e8859e108281ee7971084bfa759522	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Clostridium	 s__symbiosum
3d7262c8518a22a3315fc9fd8e2f2682	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
51ea85f428a5daaaf348c94671af386e	k__Bacteria	 p__Firmicutes	 c__Erysipelotrichi	 o__Erysipelotrichales	 f__Erysipelotrichaceae	 g__Coprobacillus	 s__
da3b4b899cecf42ff30d58718834592a	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Blautia	 s__
3f9afe17105dd0c33bf41ce08824f362	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
0dd841877cd6b20c830f4a5c7360c9ba	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
4acffcd4d2bd2707d5a65befc3d55202	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
d0097b2298c9a9d0550d27f31c049d07	k__Bacteria	 p__Cyanobacteria	 c__4C0d-2	 o__YS2	 f__	 g__	 s__
```

#### y6

```
export YY='y6'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Description1 \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv
```

Numerator has 420 OTUs

Denominator has 21

```
Feature ID	0	1	2	3	4	5	6
dd53765e55758f566920034403f46897	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Blautia	 s__
68e7e3599b8727dad8fb579a21279b55	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
eb41b1761fe5b88eb242cd3cabc8a2bb	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
869b4db9f67765ab87c9987497d4623c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales
e796a1b9569d38eef565706c91bd91fd	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
dd6c6db0c7016ecce184cc9ca510f229	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
8d5153b41e59c5ffb5a1fc080eb784cf	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
47db7fba3a4a1d17f718feb75786c18c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
77607f11477a41f99a4ef032270b77b4	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
b1e981219f9b047df26c92a9c732b845	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Faecalibacterium	 s__prausnitzii
33b04e368873d01e5d5871b9e6c79940	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
c54ecd701d322d197b9000b235740969	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
7f48094efb0f64a7ee6a2cf0c87e9636	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
7006413e51a7d9452d99564b207ea103	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
8a80b5cf4ed21ee1617aa27d0542d4a2	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
8ce212aaf7d867e84eab351d75eeecfa	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 f__Lachnospiraceae	 f__Lachnospiraceae
d56fd4c21346f6d694c7e638f6bdcb70	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
c8947944e5e63a2f031e16c4412903bd	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
4c4762e3ff3387313d4f596dd8152358	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Faecalibacterium	 s__prausnitzii
a62fa6e25e8805dafe5b628529c90023	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
95196021cd1185f319dc54ec39e5909e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
```   

#### y14

```
export YY='y14'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Description1 \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv
```
Numerator has 408 OTUs (taxa)

Denominator has 12

```
Feature ID	0	1	2	3	4	5	6
de5631cf341eed069a0d4805d5b66985	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
f6eaeb07ea32039b3d14fb82f0c8b996	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales
2f723ae2a1291cb6ed7d32ba182bd892	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
c4a48e6c0e44c439f1e6155edbcce911	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiracea e	 g__	 s__
7c6427a179cbb267c465204829fae3bf	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Clostridiaceae	 g__Clostridium	 s__
d4451ce798c2c730d419f56fa16ca60e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
2e81e4d1a54e1636af43e46fa332ab11	k__Bacteria	 p__Firmicutes	 c__Erysipelotrichi	 o__Erysipelotrichales	 f__Erysipelotrichaceae	 g__Coprobacillus	 s__
7e9bceb3c4a52554388b3d7fe500be6c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
5f5763fa626a29953839c6773446d2c4	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
e14f743e40eb38976db9d9d544548a13	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
1a069ea36017c7efc95d851374fa7382	k__Bacteria	 p__Firmicutes	 c__Bacilli	 o__Lactobacillales	 f__Lactobacillaceae	 g__Lactobacillus	 s__vaginalis
bc5991ef4c06b7bc9302a867c949570f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales
```

#### y27

```
export YY='y27'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Description1 \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv
```

Numerator has 10

```
Feature ID	0	1	2	3	4	5	6
1a069ea36017c7efc95d851374fa7382	k__Bacteria	 p__Firmicutes	 c__Bacilli	 o__Lactobacillales	 f__Lactobacillaceae	 g__Lactobacillus	 s__vaginalis
2e81e4d1a54e1636af43e46fa332ab11	k__Bacteria	 p__Firmicutes	 c__Erysipelotrichi	 o__Erysipelotrichales	 f__Erysipelotrichaceae	 g__Coprobacillus	 s__
f6eaeb07ea32039b3d14fb82f0c8b996	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales
c4a48e6c0e44c439f1e6155edbcce911	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
7c6427a179cbb267c465204829fae3bf	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Clostridiaceae	 g__Clostridium	 s__
5f5763fa626a29953839c6773446d2c4	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
2f723ae2a1291cb6ed7d32ba182bd892	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
bc5991ef4c06b7bc9302a867c949570f	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales
7e9bceb3c4a52554388b3d7fe500be6c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
d4451ce798c2c730d419f56fa16ca60e	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
```

Denominator has 2

```
Feature ID	0	1	2	3	4	5	6
e14f743e40eb38976db9d9d544548a13	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__	 s__
de5631cf341eed069a0d4805d5b66985	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
```


## OTUs that distinguish the unknown factor

#### y257

```
    export YY='y257'
    qiime gneiss balance-taxonomy \
      --i-table table-dada2.qza \
      --i-tree hierarchy.qza \
      --i-taxonomy taxonomy.qza \
      --p-taxa-level 3\
      --p-balance-name $YY \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column Unknown \
      --o-visualization ${YY}_taxa_summary.qzv

    qiime tools view ${YY}_taxa_summary.qzv
```

Numerator

```
Feature ID	0	1	2	3	4	5	6
f479c23321346918723839e33a3544f4	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
941cc85653474d3d8dfcf3fd605ab3c0	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
6d6ab7315b399f8c54b50467972c07aa	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
282cdfd2ae4301f800dd187a3ca8f559	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
e367dacc4cf58de7164775e5899d4ab2	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Ruminococcus	 s__
0b9dadade99b70f3e5f73280d3745e71	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
a230ee775999db0be47f185ac53a764c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
25b60da8f1af462e61067fc3545e8efb	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__[Ruminococcus]	 s__
07ba1b7b651f91768deae62e67afb9df	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__	 s__
350cd8e890b27c7fa7ec3a0d77dbaa9b	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
8ac806edce9a692df56a9eea816a15f0	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
6bde3249ff215fe8d72f5628f292431b	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 g__Oscillospira	 s__
0246aa933abc53fbde05942016fc2d10	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Lachnospiraceae	 g__Coprococcus	 s__
0f123b0da178612cd7dd4c9d8378ae6c	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales	 o__Clostridiales
67ed743cd7edd9fc33966deb93c8fcb2	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__Ruminococcaceae	 f__Ruminococcaceae	 f__Ruminococcaceae
```

Denominator

```
0	1	2	3	4	5	6
672e3ccf8212e58cf67f5ad92e1d4e1d	k__Bacteria	k__Bacteria	k__Bacteria	k__Bacteria	k__Bacteria	k__Bacteria	k__Bacteria  
```

The one discriminating bacteria is not well characterized.  
