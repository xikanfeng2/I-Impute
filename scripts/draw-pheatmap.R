library(pheatmap)
library(gridExtra)
library(ggpubr)

cor_data = read.csv('../experiment/heatmap/all.csv', row.names = 1)

png(filename='../experiment/heatmap/Figure3.png', width = 5000, height = 2000, res = 300)
plot_list=list()
# m = matrix(c(1:8), ncol=4)
data = read.csv('../experiment/heatmap/sim-true-counts.csv', row.names = 1)
annotation_col = read.csv('../experiment/heatmap/groups.csv', header = FALSE)
rownames(annotation_col) = colnames(data)
colnames(annotation_col)[1] = ' '
# raw data 
data = read.csv('../experiment/heatmap/sim-true-counts.csv', row.names = 1)
true = pheatmap(log10(data+1), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = '', 
        labels_row = '', annotation_col = annotation_col, annotation_legend=FALSE,
        main='Raw')
plot_list[['true']] = true[[4]]

# data = read.csv('../experiment/heatmap/counts.dca.tsv', row.names = 1)
# dca = pheatmap(log10(data+1), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = '', 
#         labels_row = '', annotation_col = annotation_col, annotation_legend=FALSE,
#         main='DCA')
# plot_list[['dca']] = dca[[4]]

data = read.csv('../experiment/heatmap/counts.saver.csv', row.names = 1)
saver = pheatmap(log10(data+1), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = '', 
        labels_row = '', annotation_col = annotation_col, annotation_legend=FALSE,
        main='SAVER')
plot_list[['saver']] = saver[[4]]

data = read.csv('../experiment/heatmap/counts.scimpute.csv', row.names = 1)
scimpute = pheatmap(log10(data+1), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = '', 
        labels_row = '', annotation_col = annotation_col, annotation_legend=FALSE,
        main='scImpute')
plot_list[['scimpute']] = scimpute[[4]]

data = read.csv('../experiment/heatmap/counts.gamma.csv', row.names = 1)
gamma = pheatmap(log10(data+1), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = '', 
        labels_row = '', annotation_col = annotation_col, annotation_legend=FALSE,
        main='C-Impute')
plot_list[['gamma']] = gamma[[4]]

data = read.csv('../experiment/heatmap/counts.gamma-saver.10.csv', row.names = 1)
iimpute = pheatmap(log10(data+1), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = '', 
        labels_row = '', annotation_col = annotation_col, annotation_legend=FALSE,
        main='I-Impute', font.main = 0.5)

plot_list[['iimpute']] = iimpute[[4]]

data = read.csv('../experiment/heatmap/sim-counts.csv', row.names = 1)
raw = pheatmap(log10(data+1), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = '', 
        labels_row = '', annotation_col = annotation_col, annotation_legend=FALSE,
        main='88.45% Dropout')
plot_list[['raw']] = raw[[4]]


ggsaver = ggscatter(cor_data, x = "SAVER", y = "Raw",
          color = "black", shape = 21, size = 0.5, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 1.5, label.y.npc = "bottom"),  
          xlab = 'Log10(count+1) in SAVER', ylab = 'Log10(count+1) in Raw',
)
plot_list[['ggsaver']] = ggsaver
ggscimput = ggscatter(cor_data, x = "scImpute", y = "Raw",
          color = "black", shape = 21, size = 0.5, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 1.5, label.y.npc = "bottom"),  
          xlab = 'Log10(count+1) in scImpute', ylab = FALSE,
)
plot_list[['ggscimput']] = ggscimput

ggcimpute = ggscatter(cor_data, x = "C.Impute", y = "Raw",
          color = "black", shape = 21, size = 0.5, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 1.5, label.y.npc = "bottom"),  
          xlab = 'Log10(count+1) in C-Impute', ylab = FALSE,
)
plot_list[['ggcimpute']] = ggcimpute

ggimpute = ggscatter(cor_data, x = "I.Impute", y = "Raw",
          color = "black", shape = 21, size = 0.5, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 1.5, label.y.npc = "bottom"),  
          xlab = 'Log10(count+1) in I-Impute', ylab = FALSE,
)
plot_list[['ggimpute']] = ggimpute

# do.call(grid.arrange,plot_list)
grid.arrange(arrangeGrob(grobs= plot_list,ncol=5))
dev.off()