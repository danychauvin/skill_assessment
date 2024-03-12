
## Taking a look at the data

For this, I will use R. I inspire myself from [[#^83e04e]]. I first need to install R and Rstudio again on my machine. See [[R]] for details about this.

Let's check the columns first:

- Forward Scatter Columns ~ Size: FSC-A, FSC-H, FSC-W (respectively Area, Height, Width)

- Side Scatter Columns ~ Granularity: SSC-A, SSC-H, SSC-W, SSC-B-A, SSC-B-H, SSC-B-W

- Fluorescence intensity channels: they all seem to be compensated in FlowJo ("FJComp"). There are 20 of them. Hence full spectrum flow cytometry.

- Time

- PCA_1 and 2, and optsne_1-2.

What about the rows: nothing. NULL.
We notice that the col names are kind of problematic. First one could get rid of all the FJComp. Second, some col names have some spaces. But let's keep it that way for now. Also have to change sample names.

Of importance: they did not substracted the autofluorescence of their cells. Because they do not autofluoresce that much (compared to macrophages for instance).

Todo
- [x] Change col names to colnames that specify the marker being investigated
- [ ] 