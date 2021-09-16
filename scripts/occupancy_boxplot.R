library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggforce)
library(grid)
library(Cairo)
args = commandArgs (trailingOnly = T)
# example line: head -5 occupancy_table/suppressed_merged_S2_to_mnase_peaks_in_open_and_closed_enhancers_spanning_lf_15_rf_15_occupancy.tsv 

# peak_877_1   Naked       17.647058823529413  closed
# peak_877_1   TF          10.294117647058822  closed
# peak_877_1   Nucleosome  72.05882352941177   closed
# peak_2105_9  Naked       85.71428571428571   open
# peak_2105_9  TF          4.081632653061225   open
figure2_theme <- function (){
    theme(plot.title=element_text( size=20 )) +
    theme(axis.title.x = element_text(colour = "black", size = 20),
          axis.title.y = element_text(colour = "black", size = 20)) +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)) +
    theme(legend.title= element_text(size = 20),
          legend.text = element_text(size = 20)) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
}


df = read.table(file("stdin"), sep = "\t", header = F, stringsAsFactors = F) 

n_open = nrow(df[df$V2 == "TF" & df$V4 == "open", ])
n_closed = nrow(df[df$V2 == "TF" & df$V4 == "closed", ])

count_closed_text = grobTree(textGrob(paste("n = ", format(n_closed, scientific = F, big.mark = ",") ), 
                             x = 0.08, y = 0.5, hjust = 0.5, rot = 90, 
                             gp = gpar(fontsize = 15)))
count_open_text = grobTree(textGrob(paste("n = ", format(n_open, scientific = F, big.mark = ",") ), 
                             x = 0.225, y = 0.5, hjust = 0.5, rot = 90, 
                             gp = gpar(fontsize = 15)))
df$V2 = factor(df$V2, levels = c("Naked", "TF", "Nucleosome"), ordered = T) 

dodge <- position_dodge(width = 1)
plt_v = ggplot(df, aes(x = V2, y = V3, fill = V4)) + geom_violin(trim = F, position = dodge) + 
        geom_boxplot(width = 0.1, position = dodge) + 
        geom_rangeframe() + theme_few() + 
        xlab("") + ylab("%Occupancy") +
        scale_fill_manual(values = alpha(c("Open" = "#386cb0", "Closed" = "#f0027f"), 0.5),
                    name = element_blank()) +
        annotation_custom(count_closed_text) + annotation_custom(count_open_text) + 
        ylim(c(0,100)) + theme (legend.position = "bottom") + 
        figure2_theme() 

Cairo::CairoPNG(args[1], width = 6, height = 4, res = 300, units = "in")
print (plt_v) 
dev.off()
        
