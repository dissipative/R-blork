# Several useful funcitons for sequence data processing
#   atgc.reduce(fasta, cutvalignment, tofile)
#   atgc.fas2nex(fas, int, gap, missing)
#   atgc.fas2phy(fas)
#   atgc.nex2fas(nex)
#   atgc.jmodeltest(fas, jmt, jre)
#   use: source("ATGC.R")
require(ape);

# Collapse sequences at given threshold.
# @param fasta - input alignment in valid FASTA format
# @param cutval - threshold for sequence deletion
# @param toFile - whether save data to file or return as funciton result
atgc.reduce <- function(fasta, cutval=0, toFile=T) {
   alignment <- read.dna(fasta, format="fasta", as.character=T);
   write(c("stay,l,excluded,l,distance"), file="cutoff.csv", append=F);
   write(c("count,a,b,distance"), file="dist.csv", append=F);
   a <- 1;
   while (a <= length(alignment[,1])) {
      b <- 1;
      while (b < length(alignment[,1])) {
         if (a == b) {
            b <- b+1;
         } else {
            di <- ( dist.dna(as.DNAbin(rbind(alignment[a,], alignment[b,])), model = "raw",
               variance = F, gamma = F,
               pairwise.deletion = T, base.freq = NULL,
               as.matrix = T) );
            cat(a,":",b, "\n",sep="");
            write( paste(a,":",b,",", row.names(alignment)[a], ",", row.names(alignment)[b], ",", di[2,1]),
               file="dist.csv",
               append=T );
            cat(di[2,1],"\n");
            if (di[2,1] == "NaN") di[2,1] <- 100;
            if (di[2,1] <= cutval) {
               length.a <- length(grep("([atgc])", as.vector(alignment[a,])));
               length.b <- length(grep("([atgc])", as.vector(alignment[b,])));
               if (length.a > length.b)
                  write(
                     paste(row.names(alignment)[a], length.a, row.names(alignment)[b], length.b, di[2,1], sep=","),
                     file="cutoff.csv",
                     append=TRUE
                     );
               if ( row.names(alignment)[a] != paste( row.names(alignment)[a],
                                                      row.names(alignment)[b],
                                                      row.names(alignment)[b],
                                                      sep=";") )
                  row.names(alignment)[a]=paste(row.names(alignment)[a],row.names(alignment)[b],sep=";");
               alignment <- alignment[-b,];
               b <- b-1;

            } else {
               write( paste(row.names(alignment)[b], length.b, row.names(alignment)[a], length.a, di[2,1], sep=","),
                  file="cutoff.csv",
                  append=T);
               if ( row.names(alignment)[b] != paste( row.names(alignment)[b],
                                                      row.names(alignment)[a],
                                                      row.names(alignment)[a],
                                                      sep=";") )
                  row.names(alignment)[b] <- paste(row.names(alignment)[b], row.names(alignment)[a], sep=";");
               alignment <- alignment[-a,];
               b <- a;
            }
         }
         b <- b+1;
      }
      a <- a+1;
   }
   if (toFile) {
      filename <- paste0( gsub(".fas","",fasta),".labeled.fas" );
      write.dna(alignment, filename, format="fasta", append=F);
   } else return(alignment);
}

# Convert alignment in FASTA to NEXUS format
# @param fasta
# @param int - interleaved whether TRUE, otherwise non-interleaved
# @param gap - gap symbol
# @param missing - the way of threating missing characters
atgc.fas2nex <- function(fasta, int=T, gap=NULL, missing=NULL) {
   input <- read.dna(fasta, format="fasta", as.character=TRUE);
   filename <- paste0(gsub(".fas", "", fasta), ".nex");
   input <- as.data.frame( t(input) );
   write.nexus.data(input, filename, format="dna", interleaved=int, gap=gap, missing=missing);
}

# Convert alignment in FASTA to non-interleaved PHYLIP format
# @param fasta
atgc.fas2phy <- function(fasta) {
   input <- read.dna(fasta, format="fasta", as.character=T);
   filename <- paste0(gsub(".fas", "", fasta), ".phy");
   input <- as.data.frame( t(input) );
   write.dna(input, filename, format="seq", indent=0, colsep="", nbcol=-1);
}

# Convert alignment in NEXUS to FASTA format
# @param nexus
atgc.nex2fas<-function(nexus) {
   input <- read.nexus.data(nexus);
   filename=paste0(gsub(".nex", "", nexus), ".fas");
   write.dna(input, filename, format="fasta");
}

# Perform the best evolution model search and prepare NEXUS file for MrBayes analysis
# @param fasta
# @param jmt - path to jModelTest jar file
# @param jre - path to java executable with '-jar' flag
atgc.jmodeltest<-function(fasta, jmt='./jModelTest.jar', jre='/usr/lib/java -jar') {

   #def params are: -s 11 -g 8 -i -f -AIC -a
   log <- system( paste0(jre," ",jmt," -d ",fasta," -s 3 -g 4 -i -f -AIC -a"), intern=T );
   log <- gsub("\t", " ",log[length(log)]);
   log <- gsub("\\s+", ",", log, perl=T);
#   names<-c("param",
#            "model",
#            "a", "c","g","t",
#            "kappa","titv",
#            "ac","ag","at","cg","ct","gt",
#            "pinv",
#            "gamma")
   params <- matrix(unlist(strsplit(log,",")), nrow=1);
   #   colnames(params)<-names

   pinv <- NULL;
   if (params[,15] == "N/A" && params[,16] == "N/A") rates<-"equal";
   if (params[,15] != "N/A") {
      pinv <- paste0("pinvarpr=fixed(",params[,15],")");
      rates <- "propinv";
   }
   if (params[,16] != "N/A") {
      shape<-paste0(" shapepr=fixed(",params[,16],")");
      rates<-"gamma";
   }
   if (params[,15] != "N/A" && params[,16] != "N/A") rates<-"invgamma";

   input <- read.dna(fasta, format="fasta", as.character=T);
   filename <- paste0(gsub(".fas","",fasta), ".nex");
   input <- as.data.frame( t(input) );
   write.nexus.data(input, filename, format="dna");
   cat(
"\n
begin mrbayes;
lset coding=alignmentl rates=",rates," nucmodel=4by4 Nst=6 nbetacat=5;
prset",pinv,shape,";
prset statefreq=dirichlet(",params[3],",",params[4],",",params[5],",",params[6],") revmatpr=dirichlet(",params[9],",",params[10],",",params[11],",",params[12],",",params[13],",",params[14],");
outgroup 1;
mcmc ngen=5000000 nruns=4 nchains=4 temp=0.200  swapfreq=1 nswaps=1 samplefreq=1000 mcmcdiagn=Yes minpartfreq=0.1 relburnin=Yes burninfrac=0.15 stoprule=No Savebrlens=Yes;
sump burnin=150 nruns=4;
sumt burnin=150 nruns=4 ntrees=1 contype=alignmentlcompat;
end;",file=filename, sep="", append=T);

cat("The best model is:",params[2]);
}