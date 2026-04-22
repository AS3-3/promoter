rm(list = ls())

library(ape)
library(readxl)
library(vegan)

# 读取数据，把基因id变成行名，并把数据转换成数值矩阵
dat <- read_excel('0413-40-矩阵.xlsx')
dat <- as.data.frame(dat)
rownames(dat) <- dat[, 1]
dat <- dat[, -1]
mat <- as.matrix(dat)

# 构建无根 NJ 树，也就是输出的树
set.seed(123)
dist_matrix <- vegdist(mat, method = "jaccard")
nj_tree <- nj(dist_matrix)

# Bootstrap 分析
B <- 1000
boot_result <- boot.phylo(
  phy = nj_tree,
  x = mat,
  FUN = function(xx) {
    d <- vegdist(xx, method = "jaccard")
    nj(d)
  },
  B = B,
  block = 1,
  trees = FALSE
)

# 外类群定根
cat("全部基因ID为：")
print(rownames(mat))

outgroup <- "IcMYB1"
if (!outgroup %in% nj_tree$tip.label) {
  stop("指定的外类群 '", outgroup, "' 不在树中，请检查名称是否正确。")
}

rooted_tree <- root(nj_tree, outgroup = outgroup, resolve.root = TRUE)

# 添加 Bootstrap 支持率
n_nodes <- Nnode(rooted_tree)
rooted_tree$node.label <- character(n_nodes)

if (length(boot_result) == Nnode(nj_tree)) {
  rooted_tree$node.label[-1] <- round(boot_result * 100 / B)
} else {
  warning("Bootstrap 结果长度与无根树内部节点数不匹配，请检查。")
}
rooted_tree$node.label[1] <- ""

# 保存 Newick 文件
write.tree(rooted_tree, file = "IbMYB1元件有根树-40-0413.nwk")

# 绘制并保存 PDF
pdf("NJ_tree_bootstrap_lengths.pdf", width = 14, height = 12)
par(mar = c(3, 3, 5, 3))

plot(rooted_tree, 
     main = "NJ tree with bootstrap values and branch lengths (Jaccard)", 
     cex = 0.8, 
     label.offset = 0.005)

# Bootstrap 支持率
nodelabels(rooted_tree$node.label, 
           adj = c(1.2, -0.2), 
           frame = "none", 
           cex = 0.7, 
           col = "blue")

# 分支长度数值
edgelabels(text = round(rooted_tree$edge.length, 3), 
           adj = c(0.5, -0.5), 
           frame = "none", 
           cex = 0.6, 
           col = "darkred")

add.scale.bar(cex = 0.8)
dev.off()

