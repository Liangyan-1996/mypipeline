#!/public/home/liunangroup/liangyan/software/miniconda3/envs/deseq2/bin/Rscript

# Rscript star_DE.R -c a.tab,b.tab,c.tab -t d.tab,e.tab,f.tab

library(argparse)

# 创建解析器对象
parser <- ArgumentParser()

# 添加参数
parser$add_argument("-l","--list",type = "character",
                    help = "sample list in csv format; one line for each condition; the first line is the control condition",
                    action = "store", default = NULL)
parser$add_argument("-s","--strandness",type = "character",
                    help = "RNA-seq library strandness. NON; FR; RF",
                    default = "NON")
# 解析命令行参数
args <- parser$parse_args()

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)


ReadCounts = function(tab,strandness){
  strand_id = c(1,2,3)
  names(strand_id) = c("NON","FR","RF")

  df = read.table(tab,sep="\t",header=F,row.names=1)
  df = df[-1:-4,] %>% select(strand_id[strandness])
  colnames(df) = gsub("\\.ReadsPerGene\\.out\\.tab", "\\1", tab)
  return(df)
}

df = lapply(c(args$control,args$treatment),function(x) {ReadCounts(x,args$strandness)})
df = do.call(cbind, lapply(df, as.matrix))

coldata = data.frame(row.names = colnames(df),
                     condition = c(rep("control",length(args$control)),rep("treatment",length(args$treatment))))

dds <- DESeqDataSetFromMatrix(countData = df, colData = coldata, design= ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition', 'treatment', 'control'))
summary(res)

res = as.data.frame(res)
res[,'sig'] <- 'none'
res[which(res$log2FoldChange >= 1 & res$padj < 0.01),'sig'] <- 'up'
res[which(res$log2FoldChange <= -1 & res$padj < 0.01),'sig'] <- 'down'

exp = read.table("1.expr_TPM.txt",sep="\t",header=T,row.names=1)
res = merge(res,exp,by = "row.names", all=F)
res <- res[order(res$padj,
                 res$log2FoldChange,
                 decreasing = c(FALSE, TRUE)),]
write.table(res,"1.DESeq2.txt",sep="\t",quote=F,row.names=F)

# BiocManager::install('EnhancedVolcano')
cc = rep("grey",nrow(res))
names(cc) = res$sig
cc["up"] = "#d73027"
cc["down"] = "#4575b4"
unique(cc)
unique(names(cc))

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = res$Gene,
                x = 'log2FoldChange',
                y = 'padj',
                ylab = bquote(~-Log[10]~ 'padj'),
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize=2,
                colCustom = cc)
ggsave("1.vocano.pdf",width=12.4,height=11.4)
