library(DiffBind)
setwd('/public/home/shidetong/projects/guolab/now_bow/difbind')
tamoxifen <- dba(sampleSheet='dif.csv')#读文件，反馈信息列表
plot(tamoxifen)#一个相关性热图
tamoxifen <- dba.count(tamoxifen)#统计reads数目和反馈FRiP

info <- dba.show(tamoxifen)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))

#数据归一化处理
tamoxifen <- dba.normalize(tamoxifen)#
norm <- dba.normalize(tamoxifen, bRetrieve=TRUE)
normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors,NormLibSize=round(norm$lib.sizes/norm$norm.factors))
#在做差异之前，先构建模型
tamoxifen <- dba.contrast(tamoxifen,reorderMeta=list(Condition="Responsive"))

#差异分析
tamoxifen <- dba.analyze(tamoxifen)
dba.show(tamoxifen, bContrasts=TRUE)
#检索差异位点
tamoxifen.DB <- dba.report(tamoxifen)
tamoxifen.DB
sum(tamoxifen.DB$Fold>0)
sum(tamoxifen.DB$Fold<0)
