import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

mrna_raw=pd.read_csv(r"C:\Users\mitteam\Desktop\GBM\Data\mrna_raw.txt", sep="\t", index_col=0)
meta_data=pd.read_csv(r"C:\Users\mitteam\Desktop\GBM\Data\meta_data.txt", sep="\t")

print(meta_data.head())
print(mrna_raw.head())

#transpose so rows=samples and cols=genes
mrna_T=mrna_raw.T
print(mrna_T.head())

tsne= TSNE(n_components=2, random_state=50, perplexity=5)
tsne_result=tsne.fit_transform(mrna_T)

tsne_df=pd.DataFrame(tsne_result, columns=["t-SNE axis 1", "t-SNE axis 2"])
tsne_df["sample_type"]=meta_data["sample_type"].values
tsne_df["sample_name"] = meta_data["sample_name"].values

#plot
plt.figure(figsize=(8,6))
Normal=tsne_df[tsne_df["sample_type"]=="Normal"]
Tumor= tsne_df[tsne_df["sample_type"]=="Tumor"]

plt.scatter(Normal["t-SNE axis 1"], Normal["t-SNE axis 2"], color="blue", label="Normal", s=100)
plt.scatter(Tumor["t-SNE axis 1"], Tumor["t-SNE axis 2"], color="red", label="Tumor", s=100)
plt.xlabel("t-SNE axis 1")
plt.ylabel("t-SNE  axis 2")
plt.title("t-SNE plot")
plt.legend()
plt.tight_layout()

plt.savefig(r"C:\Users\mitteam\Desktop\GBM\Results\my_t-SNE_plot.png", dpi=300)

plt.show()









