# data import 
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
    ref.AGI.todo <- read.table("ref.AGI.raw.txt", sep = "\t", header = T, na.strings = "")  
    ### set class
    #ref.AGI.todo[, 2] <- as.character(ref.AGI.todo[, 2])
    #ref.AGI.todo[, 3] <- as.character(ref.AGI.todo[, 3])
}

# matchup AGI between 2011MSB & Affy
matchup <- function(seq.data = seq.AGI, ref.data = ref.AGI.todo, seq.col = 1, ref.col = 2) {
    matchup.out <- data.frame(stringsAsFactors = FALSE)
    #names(matchup.out) <- c('Probe.Set.ID', 'AGI')
    for (i in c(1:dim(seq.data)[1])) {
        seq <- seq.data[i, seq.col]
        loc <- which(ref.data == seq); i; seq; loc
        if (length(loc) != 0) {
            col <- loc %/% dim(ref.data)[1]
            row <- loc %% dim(ref.data)[1]; col; row
            #matchup.row <- data.frame(as.character(ref.data[row, ref.col]), as.character(seq)); matchup.row
            matchup.row <- data.frame(ref.data[row, ref.col], seq); matchup.row
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
matchup(seq.AGI[c(1:50), ])
}
## run matchup
matchup.out <- matchup()
