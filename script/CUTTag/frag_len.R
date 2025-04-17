args = commandArgs(T)
frag_len = args[1]
out_pdf = args[2]

df = read.table(frag_len)
df$V2 = row.names(df)
colnames(df) = c("len","ID")
df = df[df$len<1000,]
library(ggplot2)
ggplot(as.data.frame(df),aes(x=len))+
  geom_density()+
  xlim(0, 1000)
ggsave(file = out_pdf, width=12, height=10,dpi=300)
