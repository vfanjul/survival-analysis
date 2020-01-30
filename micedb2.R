
doplots = F

packages = c("ggplot2","kinship2", "pedigree", "data.table")
packages = c("ggplot2","kinship2", "pedigree", "data.table", "beeswarm")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) install.packages(setdiff(packages, rownames(installed.packages())))  
for (i in packages) library(i, character.only = T)

if (.Platform$OS.type == "unix") {
  setwd("/Volumes/Victor/")
  # setwd("/Volumes/Victor-1/")
} else   setwd("S:/LAB_VA/LAB/Victor/")

route = "Mice db/"
newroute = paste0(route, Sys.Date(), "/")
dir.create(newroute, showWarnings = F)
## Generate database
webdb = c()
invisible(sapply(dir(paste0(route, "Web db/")), function (x) webdb <<- rbind.data.frame(webdb, data.frame(fread(paste0(route, "Web db/", x), skip = 1, encoding = "Latin-1", sep = ";", na.strings = c("N/A", "NA", "#N/A", ""))))))
webdb = webdb[!is.na(webdb[,1]),-c(1,18)]
names(webdb) = c("Id", "Sex", "DOB", "DOC", "Age", "Father", "Mother", "Line", "Genotype", "Status", "Substatus", "Cause", "Room", "Cage", "User", "Mark", "Project", "Comments")
webdb[webdb$Sex == "H","Sex"] = "F"
webdb$DOB = as.Date(strptime(webdb$DOB, "%d/%m/%Y"))
webdb$DOC = as.Date(strptime(webdb$DOC, "%d/%m/%Y"))
webdb = webdb[webdb$Line != "ApoEKO-CBF-b-LysMCre" & webdb$Line != "ApoEKO-CBF-b-SM22Cre",]
webdb$Genotype = gsub("LaminA/C", "Lmna", webdb$Genotype)
webdb$Genotype = gsub("  $", "", webdb$Genotype)
webdb$Genotype = gsub("^aMHC-Cre", "Lmna (G609G/LCS)  aMHC-Cre", webdb$Genotype)
webdb[webdb$Line == "C57BL/6 (H)" | webdb$Line == "C57BL/6 (CRL)", "Genotype"] = "WT"
webdb$Genotype = gsub("G609G \\(\\+/\\+", "Lmna (G609G/G609G", webdb$Genotype)
webdb$Genotype = gsub("G609G \\(\\+/wt", "Lmna (G609G/wt", webdb$Genotype)
webdb$Genotype = gsub("G609G \\(wt/wt", "Lmna (wt/wt", webdb$Genotype)
webdb$Genotype = gsub("  ", " ", webdb$Genotype)
# webdb[webdb$Genotype == "","Genotype"] = NA


webdb = webdb[order(webdb$Cage),]
webdb$CageN = table(webdb$Cage)[webdb$Cage]
# webdb$Agegroup = c(paste(seq(0,115,5), (seq(0,115,5)+4), sep = "-"), "+120")[findInterval(webdb$Age, seq(0,120,5))]
webdb$Agegroup = seq(0,120,5)[findInterval(webdb$Age, seq(0,120,5))]
webdb$DOD = webdb$DOC
webdb$Date = Sys.Date()
erratum = data.frame(fread(paste0(route,"Erratum db.txt"), na.strings = c("N/A", "NA", "#N/A", "")))
for (i in erratum$Id) webdb[webdb$Id == i, names(erratum)[!is.na(erratum[erratum$Id == i,])]] = erratum[erratum$Id == i,names(erratum)[!is.na(erratum[erratum$Id == i,])]]

## Analyze genotype
genotypes = sort(unique(webdb$Genotype))

allelelist = data.frame(t(data.frame(strsplit(gsub("\\(", "/", unique(unlist(strsplit(gsub(" ", "", genotypes), ")")))), "/"))))
allelelist$Freq = table(allelelist[,1])[allelelist[,1]]
allelelist = allelelist[allelelist$Freq > 1,]
allalleles = sort(unique(c(paste(allelelist[,1], allelelist[,2]), paste(allelelist[,1], allelelist[,3]))))
loci = data.frame("Genotype" = genotypes)
for (i in allalleles) loci[,i] = grepl(i,gsub("\\(", "", loci$Genotype)) | grepl(i,gsub("\\([[:alnum:][:punct:]]*/", "", loci$Genotype))
for (i in allalleles) loci[loci$Genotype == "WT",i] = grepl("WT", toupper(i))

webdb = merge(webdb,loci, by = "Genotype", all.x = T)
webdb = webdb[order(webdb$Id),c(names(webdb)[-1], "Genotype")]


## Status
deadmice = data.frame(table(webdb$Cause))
deadmice$FreqD = round(deadmice$Freq/nrow(webdb[webdb$Status == "Baja",])*100,1)
alivemice = data.frame(table(webdb[webdb$Status != "Baja","Status"]))
alivemice$FreqA = round(alivemice$Freq/nrow(webdb[webdb$Status != "Baja",])*100,1)
micestatus = merge(deadmice,alivemice, by = 1:2, all = T)
micestatus = merge(micestatus, data.frame(c("Alive","Dead"), c(nrow(webdb[webdb$Status != "Baja",]), nrow(webdb[webdb$Status == "Baja",]))), by = 1:2, all = T)
micestatus$Relfreq = round(micestatus$Freq/nrow(webdb)*100,1)
names(micestatus)[1] = "Status"

## Cages
micecages = merge(merge(data.frame(table(unique(webdb[!is.na(webdb$CageN) & webdb$Status != "Reproductor" & webdb$Sex == "M", c("Cage","CageN")])$CageN)),data.frame(table(unique(webdb[!is.na(webdb$CageN) & webdb$Status != "Reproductor" & webdb$Sex == "F", c("Cage","CageN")])$CageN)), all = T, by = 1), data.frame(table(unique(webdb[!is.na(webdb$CageN) & webdb$Status == "Reproductor", c("Cage","CageN")])$CageN)), all = T, by = 1)
if (nrow(webdb[webdb$Status != "Reproductor" & webdb$Sex == "M",]) == 0) micecages$Male.Stock = 0
if (nrow(webdb[webdb$Status != "Reproductor" & webdb$Sex == "F",]) == 0) micecages$Female.Stock = 0
if (nrow(webdb[webdb$Status == "Reproductor",]) == 0) micecages$Breeding = 0
names(micecages) = c("CageN", "Male.Stock", "Female.Stock", "Breeding")
micecages = merge(micecages,data.frame("Total", sum(micecages$Male.Stock, na.rm = T), sum(micecages$Female.Stock, na.rm = T), sum(micecages$Breeding, na.rm = T)), all = T, by = 1:ncol(micecages))
micecages$Total = rowSums(micecages[,2:4], na.rm = T)
micecages[is.na(micecages)] = 0

## Lines
micelines = merge(data.frame(table(webdb$Line)), data.frame(table(webdb[webdb$Status == "Baja","Line"])), all = T, by = 1)
micelines = merge(micelines, data.frame(table(webdb[webdb$Status != "Baja","Line"])), all = T, by = 1)
micelines = merge(micelines, data.frame(table(webdb[webdb$Status == "Reproductor","Line"])), all = T, by = 1)
if (nrow(webdb[webdb$Status == "Reproductor",]) == 0) micelines$BreedingM = 0
names(micelines) = c("Line", "Mice", "DeadM", "AliveM", "BreedingM")
micelines$StockM = micelines$AliveM - micelines$BreedingM
micelines = merge(micelines, data.frame(table(unique(webdb[!is.na(webdb$CageN),c("Cage","Line")])$Line)), all = T, by = 1)
micelines = merge(micelines, data.frame(table(unique(webdb[!is.na(webdb$CageN) & webdb$Status != "Reproductor",c("Cage","Line")])$Line)), all = T, by = 1)
names(micelines)[7:8] = c("Cages", "StockC")
micelines$BreedingC = micelines$Cages - micelines$StockC
micelines[is.na(micelines)] = 0

## Cage log
cagelog = data.frame(fread(paste0(route,"Mice cages log.txt"), encoding = "UTF-8", header = T, na.strings = "N/A"))
cagelog$Datecage = as.Date(cagelog$Datecage)
cagelog = rbind.data.frame(cagelog, cbind.data.frame(webdb[webdb$Status != "Baja",c("Id","Cage")], data.frame("Datecage" = c(rep(Sys.Date(), nrow(webdb[webdb$Status != "Baja",]))))))
cagelog = cagelog[order(cagelog$Datecage, decreasing = T),]
cagelog = cagelog[!duplicated(cagelog[,1:2]),]


## Population pyramids
mlines = sort(unique(webdb[,"Line"]))
ageblocks = c(paste(seq(0,115,5), (seq(0,115,5)+4), sep = "-"), "+120")

# i = mlines[1]
# i = "FAB Cre"
# i = "FABP4-Cre"
# i = "all"
if (doplots) { for (i in mlines) {
  tryCatch({
    if (!all(webdb[webdb$Line == i,"Status"] == "Baja")) {
      poppyr = data.frame(table(webdb[webdb$Line == i & webdb$Status != "Baja",c("Agegroup", "Sex", "Genotype", "Line")]))
      poppyr$Ageblock = paste(poppyr$Agegroup, (as.numeric(as.vector(poppyr$Agegroup)) + 4), sep = "-")
      poppyr[poppyr == "120-124"] = "+120"
      poppyr2 = data.frame(table(webdb[webdb$Line == i & webdb$Status != "Baja",c("Agegroup", "Sex")]))
      
      ggplot(poppyr, aes(x = Ageblock, y = Freq, fill = Genotype)) +
        geom_bar(data = poppyr[poppyr$Sex == "F" & order(poppyr$Genotype),], stat = "identity", position = "stack") +
        geom_bar(data = poppyr[poppyr$Sex == "M" & order(poppyr$Genotype),], stat = "identity", position = "stack", mapping = aes(y = -Freq)) +
        labs(x = "Age (weeks)", y = "Population", title = i) +
        coord_flip() +
        scale_x_discrete(limits = c(ageblocks, "")) +
        annotate("text", x = length(ageblocks) + 1, y = max(poppyr2$Freq)/2, label = "Female") +
        annotate("text", x = length(ageblocks) + 1, y = -max(poppyr2$Freq)/2, label = "Male") +
        theme(legend.text = element_text(size = 8)) +
        ylim(-max(poppyr2$Freq),max(poppyr2$Freq))
      ggsave(paste0(newroute, "Population pyramid ", gsub("[[:punct:]]",".", i), " .pdf"), width = 7, height = 7)
    }
  }, error = function (e) cat("ERROR :", i, " ", conditionMessage(e), "\n"))
}}

## Family trees
i = mlines[1]
for (i in c(mlines, "all")) {
  tryCatch({
    if (i == "all") {
      pre = webdb[is.na(webdb$Cause) | webdb$Cause != "Alta por error" ,c("Id","Father","Mother","Sex","Status", "Genotype", allalleles)]
    } else pre = webdb[webdb$Line == i & (is.na(webdb$Cause) | webdb$Cause != "Alta por error") ,c("Id","Father","Mother","Sex","Status", "Genotype", allalleles)]
    # pre = pre[pre$Father != "" & pre$Mother != "",]
    pre = pre[!is.na(pre$Father) & !is.na(pre$Mother),]
    parentfun = function(x) setdiff(unique(pre[!is.na(pre[,x]),x]), unique(pre[pre$Sex == ifelse(x == "Father", "M", "F"),"Id"]))
    parent = rbind.data.frame(merge(webdb[,c("Id","Father", "Mother", "Sex", "Status", "Genotype", allalleles)], parentfun("Father"), by = 1),
                              merge(webdb[,c("Id","Father", "Mother", "Sex", "Status", "Genotype", allalleles)], parentfun("Mother"), by = 1))
    if (nrow(parent) != 0) parent[,2:3] = NA
    pre = rbind.data.frame(pre, parent)
    if (length(parentfun("Father")) != 0) pre = merge(pre,data.frame(parentfun("Father"), NA, NA, "M"), by = 1:4, all = T)
    if (length(parentfun("Mother")) != 0) pre = merge(pre,data.frame(parentfun("Mother"), NA, NA, "F"), by = 1:4, all = T)
    pre$Dead = pre$Status == "Baja" | is.na(pre$Status)
    pre[grep("CX", pre$Id), "Genotype"] = "WT"
    if (any(duplicated(pre$Id))) {
      dupid = unique(pre[duplicated(pre$Id),"Id"])
      for (j in dupid) {
        pre = pre[pre$Id != j,]
        pre[(pre$Father == j | pre$Mother == j) & !is.na(pre$Father),c("Father", "Mother")] = NA
      }
    }
    if (nrow(pre) > 0 & doplots) {
      ped = with(pre, pedigree(Id,Father,Mother, Sex, rep(1,nrow(pre)), Dead))
      # if (i == "all") {
      # pdf(paste0(newroute, "Pedigree all.pdf"), 50, 10)
      # } else pdf(paste0(newroute, "Pedigree ", gsub("[[:punct:]]",".", i), ".pdf"), 50, 10)
      peddim = table(countGen(pre[order(orderPed(pre[,1:3])),1:3]))
      
      pdf(paste0(newroute, "Pedigree ", gsub("[[:punct:]]",".", i), ".pdf"), max(peddim)/2, length(peddim)*3/2)
      # pdf(paste0(newroute, "Pedigree ", gsub("[[:punct:]]",".", i), ".pdf"), max(peddim)/1.5, length(peddim)*3/1.5)
      par("usr" = c(0, max(peddim)/2, 0, length(peddim)*3/2))
      plot(ped, cex = 0.5, width = max(peddim), col = rainbow(length(sort(unique(as.numeric(factor(pre$Genotype)))))+2)[as.numeric(factor(pre$Genotype))])
      dev.off()
      legtext = sort(unique(pre$Genotype))
      if (length(legtext) > 0) {
        #   if (i == "all") {
        #     pdf(paste0(newroute, "Pedigree all legend.pdf"), 7, 7)
        #   } else pdf(paste0(newroute, "Pedigree ", gsub("[[:punct:]]",".", i), " legend.pdf"), 7, 7)
        pdf(paste0(newroute, "Pedigree ", gsub("[[:punct:]]",".", i), " legend.pdf"), max(strwidth(legtext, units = "inches")) + 2/5, (length(legtext) + 1)/5)
        par(mar = rep(0,4))
        plot(1, type = "n", axes = F, xlab = "", ylab = "")
        legend("center", pch = 16, cex = 1, legend = legtext, col = rainbow(length(sort(unique(as.numeric(factor(pre$Genotype)))))+2)[sort(unique(as.numeric(factor(pre$Genotype))))], bty = "n")
        dev.off(); dev.off()
      }
    }
  }, error = function (e) cat("ERROR :", i, " ", conditionMessage(e), "\n"))
}

pre = pre[order(orderPed(pre)),]
pre$Inbreeding = calcInbreeding(pre[,1:3])
webdb = merge(webdb, pre[,c("Id", "Inbreeding")], by = "Id", all.x = T)

fwrite(micestatus, file = paste0(newroute, "Mice status.xls"), col.names = T, sep = "\t")
fwrite(micecages, file = paste0(newroute, "Mice cages.xls"), col.names = T, sep = "\t")
fwrite(micelines, file = paste0(newroute, "Mice lines.xls"), col.names = T, sep = "\t")
for (i in c(route,newroute)) {
  fwrite(cagelog, file = paste0(i, "Mice cages log.txt"), col.names = T, sep = "\t")
  fwrite(webdb, file = paste0(i, "Mice database.xls"), col.names = T, sep = "\t")
}



