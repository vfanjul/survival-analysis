time0 = proc.time()

project = "Heart progerin" # mir29   Heart progerin   Heart lamin  ISO Challenge  PCTX Treatment  Zmpste-Rankl  HGPS Amanda  DBU Alberto
plotcensored = F # T to plot censored data in survival curves

## Import data
if (.Platform$OS.type == "unix") setwd("/Volumes/Victor/") else setwd("S:/LAB_VA/LAB/Victor/")
baseroute = paste0(project, " project/")

source("Methodology/Scripts/utility.r")

rconfigsurv = data.frame(fread(paste0(baseroute, "Design/rconfigsurv.txt"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "")))
rconfigsurv = rconfigsurv[rconfigsurv$Run == 1,-1]
invisible(sapply(names(rconfigsurv), function (x) if (class(rconfigsurv[,x]) == "integer" & x != "Maxage") rconfigsurv[,x] <<- as.logical(rconfigsurv[,x])))


mdb = data.frame(fread("Mice db/Mice database.xls", encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "", "nan")))
mdb = mdb[,c("Id", "Sex", "DOB", "Genotype", "Status", "Cause", "DOD", "Date")]
mdb$Id = toupper(mdb$Id)
mdb = mdb[!is.na(mdb$DOB),]
mdb$DOB = datefun(mdb$DOB)
mdbd = mdb[!is.na(mdb$DOD),]
mdbd$DOE = datefun(mdbd$DOD)
mdba = mdb[is.na(mdb$DOD),]
mdba$DOE = as.Date(mdba$Date)
mdb = rbind.data.frame(mdba,mdbd)

exp = 1
for (exp in 1:nrow(rconfigsurv)) {
  for(opt in 1:ncol(rconfigsurv)) assign(names(rconfigsurv)[opt], rconfigsurv[exp,opt])
  Date = as.Date(Date)
  tryCatch({
    designdata = data.frame(fread(paste0(baseroute, "Design/", study, ".txt"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "")))
    designdata$Id = toupper(designdata$Id)
    designdata$DOB = datefun(designdata$DOB)
    designdata$Start.Date = datefun(designdata$Start.Date)
    designdata$End.Date = datefun(designdata$End.Date)
    
    
    
    expdata = merge(designdata, mdb[,-4], all.x = T, sort = F)
    # expdata = expdata[expdata$DOB >= expdata$Start.Date & expdata$DOE <= expdata$End.Date,]
    # expdata = expdata[expdata$DOB <= Date,]
    # expdata$Date = format(Date, "%d/%m/%Y")
    # expdata$Date[!is.na(expdata$DOD) ] = expdata$DOD[!is.na(expdata$DOD)]
    # expdata$Date = datefun(expdata$Date)
    # expdata$Age = as.numeric(age_calc(expdata$DOB,expdata$Date, "days"))/7
    
    
    expdata$Event = ifelse(expdata$Status == "Baja", ifelse(expdata$Cause == "Encontrado muerto" | expdata$Cause == "Fin de estudio", 1, 0), 0)
    expdata$Censored = ifelse(expdata$Cause == "Encontrado muerto" | expdata$Cause == "Fin de estudio", 0, 1)
    expdata = expdata[!is.na(expdata$DOE),]
    expdata$Age = as.numeric(age_calc(expdata$DOB,expdata$DOE, "days"))/7
    
    
    for (i in sort(unique(expdata$Condition.Id))) expdata$Condition.Id[expdata$Condition.Id == i] = which(sort(unique(expdata$Condition.Id)) == i)
    for (i in sort(unique(expdata$Condition))) expdata[,i] = expdata$Condition == i
    
    dir.create(paste0(baseroute, "Results"), showWarnings = F)
    dir.create(paste0(baseroute, "Results/", study), showWarnings = F)
    route = paste0(baseroute, "Results/", study, "/Survival/")
    dir.create(route, showWarnings = F)
    
    
    ## Stats
    groups = unique(expdata$Condition[order(expdata$Condition.Id)])
    factors = c()
    if (sex_factor) factors = c(factors, "Sex")
    
    survstats = data.frame(matrix(ncol = 1))[-1,]
    for (i in 1:length(groups)) {
      cox = summary(coxph(as.formula(paste("Surv(Age, Event) ~", paste(c(factors, groups[-i]), collapse = " + "))), expdata))
      assign(paste0("survstats", i), cbind.data.frame("Comparison" = c(row.names(cox[[7]])[1:(length(row.names(cox[[7]])) - length(groups[-i]))], paste0(groups[-i], sep = paste0("vs", groups[i]))),
                                                      cox[[7]], stringsAsFactors = F))
      survstats = rbind.data.frame(survstats,get(paste0("survstats",i)), make.row.names = F, stringsAsFactors = F)
    }
    survstats[,-1] = round(survstats[,-1],5)
    survstats = survstats[!duplicated(abs(survstats[,c(2,4:6)])),]
    if (any(factors == "Sex")) survstats[survstats$Comparison == "SexM","Comparison"] = "MalevsFemale"
    if (any(factors == "Seasons")) survstats[survstats$Comparison == "Seasons","Comparison"] = "WintSprivsSumFall"
    names(survstats)[6] = "p.value"
    survstats$Sig = sapply(survstats[,6], function (x) ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", round(x, 4)))))
    survstats$Sign = ifelse(survstats[,6] < 0.05, - sign(survstats[,2]), NA)
    
    # expdata$Season = relevel(expdata$Season, "Spring")
    # expdata$Season = relevel(expdata$Season, "Winter")
    # expdata$Season = relevel(expdata$Season, "Summer")
    
    
    ## Stats for plots
    colors = colfun(length(groups))
    if ("Color" %in% names(expdata)) colors = unique(expdata$Color[order(expdata$Condition.Id)])
    
    plotstats = survstats[,c("Comparison","p.value","Sig")]
    levelsfactors = c()
    for (i in factors) levelsfactors = c(levelsfactors, levels(as.factor(expdata[,i])))
    if (length(factors) > 0) plotstats = plotstats[-1:-(length(levelsfactors)-length(factors)),]
    plotstats$G1 = sapply(gsub("vs[[:alnum:][:punct:]]*", "", plotstats$Comparison), function(x) which(groups == x))
    plotstats$G2 = sapply(gsub("[[:alnum:][:punct:]]*vs", "", plotstats$Comparison), function(x) which(groups == x))
    plotstats$P1 = sapply(1:nrow(plotstats), function (x) min(plotstats$G1[x],plotstats$G2[x]))
    plotstats$P2 = sapply(1:nrow(plotstats), function (x) max(plotstats$G1[x],plotstats$G2[x]))
    plotstats$Mean = rowMeans(data.frame(plotstats$P1, plotstats$P2))
    plotstats$Length = plotstats$P2 - plotstats$P1
    plotstats = merge(plotstats, plotsigaxes[plotsigaxes$Ngroups == length(groups),], by = c("P1","P2"), all.x = T)
    heightsigaxes = data.frame("Height" = sort(unique(plotstats[plotstats$p.value < 0.05, "Height"])))
    if (nrow(heightsigaxes) > 0) {
      heightsigaxes$Heightsig = seq(1, length(sort(unique(plotstats[plotstats$p.value < 0.05, "Height"]))), 1) - 1
      plotstats = merge(plotstats, heightsigaxes, by = c("Height"), all.x = T)
    } else plotstats$Heightsig = plotstats$Height  
    if (onlyplotsigcomps) {
      plotaxes = plotstats[!is.na(plotstats$Heightsig),]
    } else {
      plotaxes = plotstats[plotstats$p.value < 0.05 & !is.na(plotstats$Heightsig),]
      plotaxes$Height = plotaxes$Heightsig
    }
    
    ## Export freqs &  stats table and plot
    freqdata = table(expdata$Condition)
    if (length(factors) > 0) for (i in factors) freqdata = cbind(freqdata, table(expdata[, c("Condition", i)]))
    freqdata = data.frame(names(table(expdata$Condition)), freqdata)
    names(freqdata)[1:2] = c("Condition", "Total")
    
    fwrite(freqdata, file = paste0(route, "frequencies.txt"), row.names = F, sep = "\t")
    fwrite(survstats, file = paste0(route, "survstats.xls"), row.names = F, sep = "\t")
    fwrite(expdata, file = paste0(route, "expdata.xls"), row.names = F, sep = "\t")
    
    pdf(paste0(route, "Condition stats.pdf"), (length(groups) + 1)/5, (max(plotaxes[,"Height"]) + 3)/5)
    par(bty = "n", cex.axis = 0.7, mar = c(0,0,max(plotaxes[,"Height"]) + 1,0))
    plot(1:length(groups), rep(0, length(groups)), pch = 16, cex = 1, col = colors, xaxt = "n", yaxt = "n", xlab ="", ylab = "", ylim = c(0,0))
    if (nrow(plotaxes) > 0) for (j in 1:nrow(plotaxes)) {
      axis(side = 3, line = plotaxes[j,"Height"], at = plotaxes[j,c("P1","P2")] + c(0.1,-0.1), tck = 0, labels = F, lwd = 2)
      mtext(side = 3, line = plotaxes[j,"Height"], at = plotaxes[j,"Mean"], plotaxes[j,"Sig"], cex = 0.8)
    }
    dev.off()  
    

    ### Survival plots
    i = "Condition.Id"
    for (i in c(factors, "Condition.Id")) {
      pdf(paste0(route, i, " survial plot.pdf"), 10/0.9, 4/0.9, pointsize = 24, useDingbats = F)
      par(bty = "l", cex.axis = 0.75, mar = c(3,3,1,1), mgp = c(1.5,0.25,0), tck = - 0.03)
      plotlim = Maxage
      if (max(expdata$Age) < Maxage) plotlim = max(expdata$Age)
      # plot(c(0, max(expdata$Age)), c(0,1), xlab = "Age (weeks)", ylab = "Survival Probability", xaxs = "i", yaxs = "i", xaxp = c(0, round(max(expdata$Age)/10,0)*10, round(max(expdata$Age)/10,0)*5), col = 0)
      plot(c(0, max(expdata$Age)), c(0,1), xlab = "Age (weeks)", ylab = "Survival Probability", xaxs = "i", yaxs = "i", xaxp = c(0, round(max(expdata$Age)/10,0)*10, round(max(expdata$Age)/10,0)*5), xlim = c(0, plotlim), col = 0)
      abline(h = 0.5, col = 1, lty = 3)
      for (j in 1:length(sort(unique(expdata[, i])))) lines(survfit(Surv(Age, Event) ~ 1, expdata[as.numeric(as.factor(expdata[,i])) == j,]), mark.time = plotcensored, conf.int = F, col = colors[j], lwd = 3)
      dev.off()
      if (i == "Condition.Id") {
        legtext = groups
      } else legtext = (sort(unique(expdata[, i])))
      pdf(paste0(route, i, " legend.pdf"), max(strwidth(legtext, units = "inches")) + 2/5, (length(legtext) + 1)/5)
      par(mar = rep(0,4))
      plot(1, type = "n", axes = F, xlab = "", ylab = "")
      legend("center", pch = 16, cex = 1, legend = legtext, col = colors, bty = "n")
      dev.off()  
    }
  }, error = function (e) cat("ERROR :", study, " ", technique, " ", conditionMessage(e), "\n"))
}