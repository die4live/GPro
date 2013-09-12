# data import
import.raw <- read.table(file.choose())
## AGI of interest from 2011MSB -> seq.AGI
{
    seq.AGI <- read.table("seq.AGI.txt", sep = "\t", header = T)
    ### set class
    seq.AGI$AGI <- as.character(seq.AGI$AGI)
    seq.AGI$GPA1.Prey <- as.factor(seq.AGI$GPA1.Prey)
    seq.AGI$AGB1.Prey <- as.factor(seq.AGI$AGB1.Prey)
}
## AGI ref from Affy -> ref.AGI
{
    ref.AGI <- read.table("ref.AGI.raw.txt", sep = "\t", header = T, na.strings = "")  
    ### set class
    #ref.AGI[, 2] <- as.character(ref.AGI[, 2])
    #ref.AGI[, 3] <- as.character(ref.AGI[, 3])
}
## expr data from 2008PM -> expr.raw && expr.signal
{
    expr.raw <- read.table("expr.raw.txt", sep = "\t", header = T)
    expr.signal <- expr.raw[, c(1,2,3,6,9,12,15,18,21,24)]  ## extract signal data
}

# FUNC matchup AGI between 2011MSB & Affy
> FUNC matchup() aims to matchup cells in columns (search from <seq.col> in <ref.col>) between two dataframe (<seq.data> and <ref.data>) then outputs rows matchuped containing interested columns (<extra.col>)
matchup <- function(seq.data = seq.AGI, ref.data = ref.AGI, seq.col = 1, ref.col = 2, extra.col = 1) {
    matchup.out <- data.frame(stringsAsFactors = FALSE)
    #names(matchup.out) <- c('Probe.Set.ID', 'AGI')
    for (i in c(1:dim(seq.data)[1])) {
        seq <- seq.data[i, seq.col]
        loc <- which(ref.data == seq); i; seq; loc
        if (length(loc) != 0) {
            col <- loc %/% dim(ref.data)[1]
            row <- loc %% dim(ref.data)[1]; col; row
            #matchup.row <- data.frame(as.character(ref.data[row, ref.col]), as.character(seq)); matchup.row
            matchup.row <- data.frame(ref.data[row, ref.col], seq, ref.data[row, extra.col]); matchup.row
            matchup.out <- rbind(matchup.out, matchup.row); matchup.out
        }
    }
    return(matchup.out)
}
## test matchup
{
    seq.test <- data.frame(c("ATMG00880", "AT2G07768", "ATMG00960"), c(1,0,1))  ## ref.data[5,3],[6,3][6,4]
    seq.test[,1] <- as.character(seq.test[,1])
    matchup(seq.test)
    rm(seq.test)
    matchup(seq.AGI[c(1:50), ])
}
## run matchup -> matchup.out
{
    matchup.out <- matchup()
    names(matchup.out) <- c('Affy', 'AGI', 'Set.Num')
    #matchup.out_id <- matchup(seq.data = matchup.out, ref.col = 1)  ## add Set.Num !fail
}

# expr comp signal
## matchup expr gpro of matchup.out to 2008PM data -> expr.out
expr.out <- matchup(seq.data = matchup.out, ref.data = expr.signal, seq.col = 3, ref.col = 1, extra.col = 2:10)
## FUNC to clean duplicate Set.Num # errors caused unknown actually
> FUNC clean.dup aims to cleanup duplicate rows containing same cells bewteen two columns (<ref.col> and <dup.col>)
clean.dup <- function(data = expr.out, ref.col = 1, dup.col = 2)
{
    expr = data.frame()
    for (i in c(1:dim(data)[1])) {
        if (expr.out[i, ref.col] == expr.out[i, dup.col]) {
            expr <- rbind(expr, expr.out[i,])
        }
    }
    return(expr)
}
## run clean.dup -> expr && geneID
{
    expr <- clean.dup()
    expr <- expr[, 2:dim(expr)[2]]
    rownames(expr) <- seq(1:dim(expr)[1])
    names(expr) <- c('Set.Num', 'Affy', 'GCI', 'GAI', 'MCI', 'MAI', 'GCN', 'GAN', 'MCN', 'MAN')
    expr <- expr[, c(1, 3:dim(expr)[2])]
}
## Notations of expr treatments
{
    G :: guard cell
    M :: mesophyll cell
    C :: control treatment
    A :: ABA treatment  *100 uM*
        I :: with inhibitors *Actinomycin and Cordycepin were added during protoplast isolation.*
        N :: no inhibitor
}
## expr.out export -> expr.out && expr.out.txt
{
    expr.out <- cbind(matchup.out,expr)
    expr.out <- expr.out[,-4]
    write.table(expr.out, "./expr.out.txt", sep="\t")
}

# expr signal explore
## order by Set.Num
var <- expr[order(expr$Set.Num), ]
`x[order(x[,1],x[,2]),]` for multi-col order
row.names(var) <- NULL
## calculate signal values between two columns as percentage
> FUNC comp() compares values of two column names (<comp.col>, <target.col>) output a column (<name>)
comp <- function(comp.col, target.col, name = 'comp', data = var, rounds = 3)
{
    c.col <- which(colnames(data)==comp.col)
    t.col <- which(colnames(data)==target.col)
    data <- cbind(data, round((data[, c.col] - data[, t.col]) / data[, t.col], rounds))
    colnames(data)[dim(data)[2]] <- name
    return(data)
}
> similar to comp() but output results only as a vector object
comp.col <- function(comp.col, target.col, data = var, rounds = 3)
{
    c.col <- which(colnames(data)==comp.col)
    t.col <- which(colnames(data)==target.col)
    col <- round((data[, c.col] - data[, t.col]) / data[, t.col], rounds)
    return(col)
}
### test comp()
comp('GAI', 'GCI','GvI')
## calc var between treatments
{
var <- comp('GAI', 'GCI', 'GvI')
var <- comp('MAI', 'MCI', 'MvI')
var <- comp('GAN', 'GCN', 'GvN')
var <- comp('MAN', 'MCN', 'MvN')
var <- comp('MCN', 'GCN', 'vCN')
var <- comp('MAN', 'GAN', 'vAN')
}
## explore data
{
boxplot(var[, c(10:13)])
boxplot(var[, c(14:15)])
barplot(var$GCN, names.arg = var$Set.Num)
barplot(t(as.matrix(var[, c(2:9)])), beside = T, names.arg=var$Set.Num)
barplot(t(as.matrix(var[, c(10:15)])), beside = T, names.arg=var$Set.Num)
}

# filter data
> FUNC outline filters lower and upper outliners in <data> according to column <name> and output a dataframe with <ID.name> and values
outline <- function(name, data = var, ID.col = 1, ID.name = 'Set.Num')
{
    col.n <- which(colnames(data) == name)
    col <- data[, col.n]
    margin <- (quantile(col)[4] - quantile(col)[2]) * 1.5
    lower <- quantile(col)[2] - margin
    upper <- quantile(col)[4] + margin
    names(margin) <- NULL
    names(lower) <- NULL
    names(upper) <- NULL
    out.lower <- data[which(col < lower), col.n] ; out.lower
    out.upper <- data[which(col > upper), col.n] ; out.upper
    num.lower <- numeric()
    m <- 1
    for (i in out.lower) {
        num.lower[m] <- data[which(data[, col.n] == i), ID.col]
        m <- m + 1
    }
    num.upper <- numeric()
    n <- 1
    for (i in out.upper) {
        num.upper[n] <- data[which(data[, col.n] == i), ID.col]
        n <- n + 1
    }
    lower <- rbind(num.lower, out.lower)
    lower <- rbind(lower, c(rep(paste(name, 'lower'), dim(lower)[2])))
    upper <- rbind(num.upper, out.upper)
    upper <- rbind(upper, c(rep(paste(name, 'upper'), dim(upper)[2])))
    colnames(lower) <- c(rep(paste(name, 'lower'), dim(lower)[2]))
    colnames(upper) <- c(rep(paste(name, 'upper'), dim(upper)[2]))
    rownames(lower) <- c(ID.name, 'var.value', 'var.type')
    rownames(upper) <- c(ID.name, 'var.value', 'var.type')
    out <- cbind(lower, upper)
    out <- as.data.frame(out)
    return(out)
}
## output outliners
out <- cbind(outline('GvI'), outline('MvI'), outline('GvN'), outline('MvN'), outline('vCN'), outline('vAN'))
colnames(out) <- c(1:dim(out)[2])
## matchup results to Affy or AGI
ref <- matchup.out[order(matchup.out$Set.Num), ]
rownames(ref) <- NULL
out.t <- t(as.matrix(out))
rm(expr, expr.out, expr.raw, matchup.out, out, ref.AGI, seq.AGI)
out.gene <- matchup(seq.data = out.t, ref.data = ref, seq.col = 1, ref.col = 3, extra.col = 1:2)
out.gene <- cbind(out.gene[, c(3:4)], out.t)
out.gene <- out.gene[order(out.gene$Set.Num), ]
rownames(out.gene) <- NULL
## export results & plots
write.table(out.gene, "./archive/out.gene.txt", sep="\t")