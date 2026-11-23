# MorphGate

和 GHIST 的对照写在 [docs/article.md](docs/article.md)。协议来自 Wang 等人，Nature Communications，2025。

H&E 上看不见的基因，就不拿图像去预测。MorphGate 在 HER2 阳性乳腺癌上做留一病人：先在训练集里按字母留出一位病人，只保留 Pearson 相关至少 0.10 的基因，再用斑点及其六个近邻的苏木精/伊红特征做岭回归。没通过的基因回填训练集均值，这条常数预测的相关记为 0。

测试病人不是 A 时，门控病人都是 A。流水线没有随机划分，方差排序、门控病人和岭回归都是确定的。

数据是 [Andersson 等，Zenodo 4751624](https://doi.org/10.5281/zenodo.4751624)，8 位病人、36 张切片。表达是 log1p(CP10k)。每个折先在其余病人上取方差最高的 250 个基因。

留一病人的结果：平均留下 72.1 个基因，这些基因的平均 Pearson 相关是 0.139，八个折从 0.047（病人 F）到 0.211（病人 C）。被门控丢掉的基因相关记为 0。250 个基因合在一起，平均相关是 0.037。反复留下的是胶原、HER2 通路和上皮角蛋白，例如 FN1、COL1A1、ERBB2、FASN、KRT18。

```bash
python -m morphgate --root data/her2st --out results
```

## 现有方法

综述 *Predicting spatial transcriptomics from histology images: Progress, challenges, and directions*（Exploration，2026）把四十多个模型收成两类输入，再按学习方式分开。报告的 Pearson 相关大多在 0 到 0.6，外部验证很少，训练集偏乳腺癌 Visium。标记基因和空间可变基因比全转录组容易。切片内部或连续切片的留出，会把相关抬高。

**单斑块回归。** ST-Net、DeepSpaCE、BrST-Net、HE2Gene、HistoSPACE、HEST-1k 的岭回归、STimage、Path2Space、sCellST、NSL。每个斑点单独从图像回归表达。优点是省内存，小模型往往不比大图网络差。ST-Net 的外部平均相关约 0.33，DeepSpaCE 约 0.13，HEST-1k 在 Visium 上约 0.23–0.39。缺点是看不见邻域，小队列上容易记住染色批次。

**带空间上下文的网络。** HisToGene、Hist2ST、THItoGene、TCGN、TRIPLEX、SEPAL、EGN、EGGN、ErwaNet、STFormer、HGGEP、STPath、PH2ST。用图网络或 Transformer 把邻近斑点或整张切片连起来。优点是能表示更大的组织结构。综述引用的比较里，这种上下文并没有稳定超过简单卷积网络。TRIPLEX 的外部相关约 0.12–0.14，STFormer 的外部中位数约 0.15。图网络还容易把邻域抹平，整张切片一次读入内存。

**对比学习和检索。** BLEEP、HECLIP、OmiCLIP、CMRCNet、MclSTExp、STMCL、RankByGene、UMPIRE。把图像和表达嵌到同一空间，预测时取参考斑点的表达。优点是可以找回一个真实的表达谱，而不是回归出一个平均数。缺点是查询图像一旦离开参考库，就只能抄近邻；批次效应会主导嵌入。BLEEP 在连续切片上约 0.17–0.22，OmiCLIP 在 ST-bank 上中位数约 0.39，UMPIRE 约 0.29。

**病理基础模型加轻量头。** HEST-1k、STimage、Path2Space、DeepSpot、STFlow、STPath、PH2ST、Stem、RankByGene 使用 UNI、CONCH、GigaPath、Phikon、CTransPath、Virchow、H-Optimus。优点是形态特征不用从头学，往往只要训练一个线性层。缺点是权重大、使用受限，而且这些模型是在较大图块上预训练的，换到单细胞图块时分布已经变了。几个基础模型拼在一起，增益很小。

**生成模型。** XFuse、Stem、STFlow。可以做超分辨率或按图像采样表达。计算更重。XFuse 在连续切片上报告过很高的相关，这种划分本身会抬高数字，细结构也很少有正交实验核对。

**改预测对象。** GHIST 预测细胞类型再回推表达，DeepPT、DeepPathway 一类工作预测签名或通路。形态学看得见的是细胞组成和组织区域。通路和区域比单基因稳，但不再是完整转录组。

MorphGate 取单斑块这一支，只加六个近邻的平均染色，并显式丢掉在另一位病人上没有关联的基因。它不使用病理基础模型，因此不依赖需要审批的权重。0.139 只覆盖通过门控的基因；把丢掉的基因算进去，相关是 0.037。

## 文献里的数据

反复被用作训练和外部测试的公开集放在服务器 `/mnt/storage/bzhang/morphgate/data`：

- HER2ST，36 张切片，8 个供体，约 1.2 GB。Zenodo [10.5281/zenodo.4751624](https://doi.org/10.5281/zenodo.4751624)。
- 10x Visium 人乳腺癌 Block A Section 1，约 47 MB。这是过滤后的表达矩阵，加上空间坐标和高分辨率预览图，综述和 ST-Net 用它做外部切片。
- 皮肤鳞癌处理包，约 2.5 GB，Hist2ST 等使用的 cSCC 队列。Zenodo [13117873](https://zenodo.org/records/13117873)。

HEST-1k、STimage-1K4M、ST-bank、ViSTomics-4M 各有约一千到一千四百张切片、两百万到四百万个图块。磁盘上另有正在写入的队列，这四份全集没有整库镜像。
