# rnaseqVis

[rnaseqVis app on posit Connect Cloud](https://019658d9-521a-3d9f-fd3f-8d5cd9302c19.share.connect.posit.cloud)

## Instructions

The app is designed to display gene expression data produced by RNA-seq.

### Loading data

The input page is used to load the RNA-seq data files. Data can be loaded using the buttons at the top of the side bar control panel. Users need to supply a sample file along with a count file.

#### Count File
The count file should be a tab-separated file containing genes in rows and samples in columns.
The 'GeneID' column is required. If you are using IDs that don't obviously correspond to a gene name, simply duplicate the IDs column. The other required columns are the counts for each different sample either as already normalised counts or raw counts.
If normalised count columns are provided these will be used directly (The column names must match the sample IDs with ' normalised count' appended). Otherwise count columns will be loaded into DESeq2 and normalised for library size (In this case the column names must match the sample IDs with ' count' appended). The normalised counts are then displayed. If both normalised and unnormalised counts are provided, the normalised counts will be used.

e.g.
<table class="table-bordered-centered">
  <tr>
    <th>GeneID</th>
    <th>Gene name</th>
    <th>012525 count</th>
    <th>012526 count</th>
    <th>012536 count</th>
    <th>012537 count</th>
    <th>012525 normalised count</th>
    <th>012526 normalised count</th>
    <th>012536 normalised count</th>
    <th>012537 normalised count</th>
  </tr>
  <tr>
    <td>ENSDARG00000000212</td>
    <td>krt97</td>
    <td>10</td>
    <td>12</td>
    <td>35</td>
    <td>42</td>
    <td>12.3</td>
    <td>14.6</td>
    <td>35.6</td>
    <td>45.1</td>
  </tr>
  <tr>
    <td>ENSDARG00000000567</td>
    <td>znf281a</td>
    <td>345</td>
    <td>333</td>
    <td>357</td>
    <td>365</td>
    <td>322.5</td>
    <td>343.7</td>
    <td>363.2</td>
    <td>380.0</td>
  </tr>
</table>

#### Sample File

The sample file is a tab-separated file containg at least two columns. The required columns are the sample IDs that match the (normalised) count column names in the count file. This must be the first column. The other required column should be labelled "condition" and details how the samples are divided into groups (e.g. Control and Treated). If the file contains a column labelled "sampleName" these names will be displayed on the heatmap instead of the IDs from the first column. The count columns will be reordered based on the sample file.

e.g.
<table class="table-bordered-centered">
  <tr>
    <th>sample</th>
    <th>condition</th>
    <th>sampleName</th>
    <th>batch</th>
  </tr>
  <tr>
    <td>012536</td>
    <td>Control</td>
    <td>Ctrl1</td>
    <td>A</td>
  </tr>
  <tr>
    <td>012537</td>
    <td>Control</td>
    <td>Ctrl2</td>
    <td>B</td>
  </tr>
  <tr>
    <td>012525</td>
    <td>Treated</td>
    <td>Trt1</td>
    <td>A</td>
  </tr>
  <tr>
    <td>012526</td>
    <td>Treated</td>
    <td>Trt2</td>
    <td>B</td>
  </tr>
</table>

#### Subset by Gene id

A text file of Ensembl gene ids to subset the heatmap to can be uploaded using this button. Any ids that can't be matched will be listed in a warning alert above the heatmap.

#### Transform counts

The counts can be transformed as follows

* Raw - The untransformed normalised counts (Default)
* Max Scaled - Each row is scaled to the maximum value for that row.
* log10
* Mean Centred and Scaled - Each row is centred by subtracting the mean and scaled by dividing by the standard deviation for that row.

#### Clustering

Both the rows (genes) and columns (samples) of the heatmap can be clustered. Currently, the clustering is done by hierarchical clustering of the Pearson correlation coefficients between genes/samples. Only the genes/samples currently displayed in the heatmap are used for the clustering
