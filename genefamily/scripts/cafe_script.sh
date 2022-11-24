#! ~/software/CAFE/release/cafe
date
load -i cafeinput.filtered -p 0.01 -t 10 -l log.txt
# -t, 使用的线程数
# -p, 显著性判断阈值
#the phylogenetic tree structure with branch lengths
# 这棵树可以直接用mcmctree输出的树
tree (Nnu:124.5781,((Ath:24.2897,Tha:24.2897):94.0914,(((((Bse:13.0498,Bgy:13.0498):2.2495,Bcy:15.2993):22.3403,Rap:37.6396):15.1566,Clo:52.7963):40.9378,((Mes:67.4141,Rco:67.4141):19.6264,(Spu:22.7468,(Ptr:10.8282,Peu:10.8282):11.9186):64.2937):6.6936):24.6470):6.1971)
#search for 2 parameter model
# 也可以对不同的clade设置不同的lambda值，最终的结果好像差别不是很大
lambda -s -t (1,((1,1)1,(((((1,1)1,1)1,1)1,1)1,((1,1)1,(1,(1,1)1)1)1)1)1)
report reportout1
date
