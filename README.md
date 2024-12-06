# corHex

corHex is an R package for generating correlation hexagon plots. It helps you quickly visualize correlations between multiple groups within the R environment, particularly suited for large-scale datasets.

# Installation

```
# Install devtools
install.packages("devtools")

# Install corHex package from GitHub
devtools::install_github("Chuanping-Zhao/corHex")
```
# Example Data

dtplot:
![image](https://github.com/user-attachments/assets/941e76b6-3e63-41a7-99bc-8e84190e28b2)

# Using the cor_hex Function

```
corHex::cor_hex(
dt=dtplot,
id.col="protein",
cor.method=c("pearson", "kendall", "spearman")[1],
savefile="outputfile",
singleplotsize=c(3,2.5),#width height
facetplotsize=c(3*3,2.5*3),#width height
                bin=50
)
```
