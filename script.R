# data import
import.raw <- read.table(file.choose())
## AGI of interest from 2011MSB
{
    seq.AGI <- read.table("seq.AGI.txt", sep = "\t", header = T)
    ### set class
    seq.AGI$AGI <- as.character(seq.AGI$AGI)
    seq.AGI$GPA1.Prey <- as.factor(seq.AGI$GPA1.Prey)
    seq.AGI$AGB1.Prey <- as.factor(seq.AGI$AGB1.Prey)
}
## AGI ref from Affy
{
    ref.AGI <- read.table("ref.AGI.raw.txt", sep = "\t", header = T, na.strings = "")  
    ### set class
    #ref.AGI[, 2] <- as.character(ref.AGI[, 2])
    #ref.AGI[, 3] <- as.character(ref.AGI[, 3])
}
## expr data from 2008PM
{
    expr.raw <- read.table("expr.raw.txt", sep = "\t", header = T)
    expr <- expr.raw[, c(1,2,3,6,9,12,15,18,21,24)]  ## extract signal data
}

# matchup AGI between 2011MSB & Affy
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
## run matchup
{
matchup.out <- matchup()
names(matchup.out) <- c('Affy', 'AGI', 'Set.Num')
#matchup.out_id <- matchup(seq.data = matchup.out, ref.col = 1)  ## add Set.Num !fail
}

# expr comp signal
