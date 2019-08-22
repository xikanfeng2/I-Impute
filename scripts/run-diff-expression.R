library(edgeR)
data = read.csv('GSE74672.csv', row.names = 1)
sg <- character(length(colnames(data)))
for (i in 1:length(colnames(data))){
    sg[i] <-  strsplit(colnames(data)[i], '\\.')[[1]]
}
names(sg) <- colnames(data)
sg <- factor(sg)
y <- DGEList(counts = data, group = sg)
keep <- rowSums(cpm(y)>1) >= 3;
y <- y[keep,]
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
top <- topTags(et,n=2000)
