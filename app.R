library(shiny)
library(vote)
library(shinythemes)
library(scales)
library(ggplot2)
library(MCMCpack)

#will have to think about these parameters
sample_interval = 20
burn = 1000
numcandtypes <- 2 
globalcandtypes <- list(preferred_1 = TRUE, preferred_2 = FALSE)

valid_prefix <- function(v){
    if (all(!is.na(v))){
        return(v)
    }
    if (any(!is.na(v))){
        return(v[1:min(which(is.na(v)))-1])
    }
    else{
        return(c())
    }
}
get_positions <- function(df, row, num_positions){
    list(
        minpositions = unlist(valid_prefix(df[row,3:(num_positions+2)])),
        majpositions = unlist(valid_prefix(df[row,(num_positions+3):ncol(df)]))
    )
}

multiGroupBallotsCambridgeSampler <- function(
    minority_crossover, majority_crossover,
    minority_candidates, majority_candidates,
    minority_vary_minority, minority_vary_majority,
    majority_vary_minority, majority_vary_majority,
    minority_share, numballots, ballot_type_frames){
        ballots = matrix(nrow=0,ncol=length(minority_candidates)+length(majority_candidates))
        k = length(majority_candidates)+length(minority_candidates)
        colnames(ballots) <- c(1:k)
        #maj block
        theseballots <- sapply(
            seq(floor(numballots*(1-minority_share)*(1-majority_crossover))),
            function(x){
                if (majority_vary_majority) {
                    majorder <- sample.int(length(majority_candidates))
                }
                else {
                    majorder <- c(1:length(majority_candidates))   
                }
                if (majority_vary_minority) {
                    minorder <- sample.int(length(minority_candidates))
                }
                else {
                    minorder <- c(1:length(minority_candidates))   
                }
                chosen_type <- sample(1:nrow(ballot_type_frames$Wfirst), prob=ballot_type_frames$Wfirst$p, size=1)
                positions <- get_positions(ballot_type_frames$Wfirst,chosen_type,20)
                majpositions <- positions$majpositions[1:min(length(majority_candidates),length(positions$majpositions))]
                minpositions <- positions$minpositions[1:min(length(minority_candidates),length(positions$minpositions))]
                b <- vector()
                length(b) <- k
                b[majpositions] <- majority_candidates[majorder[1:length(majpositions)]]
                b[minpositions] <- minority_candidates[minorder[1:length(minpositions)]]
                ballot <- valid_prefix(b)
                length(ballot) <- k
                ballot
            }
        )
        theseballots <- t(theseballots)
        ballots <- rbind(ballots,theseballots)
        #maj cross
        theseballots <- sapply(
            seq(floor(numballots*(1-minority_share)*majority_crossover)),
            function(x){
                if (majority_vary_majority) {
                    majorder <- sample.int(length(majority_candidates))
                }
                else {
                    majorder <- c(1:length(majority_candidates))   
                }
                if (majority_vary_minority) {
                    minorder <- sample.int(length(minority_candidates))
                }
                else {
                    minorder <- c(1:length(minority_candidates))   
                }
                chosen_type <- sample(1:nrow(ballot_type_frames$Cfirst), prob=ballot_type_frames$Cfirst$p, size=1)
                positions <- get_positions(ballot_type_frames$Cfirst,chosen_type,20)
                majpositions <- positions$majpositions[1:min(length(majority_candidates),length(positions$majpositions))]
                minpositions <- positions$minpositions[1:min(length(minority_candidates),length(positions$minpositions))]
                b <- vector()
                length(b) <- k
                b[majpositions] <- majority_candidates[majorder[1:length(majpositions)]]
                b[minpositions] <- minority_candidates[minorder[1:length(minpositions)]]
                ballot <- valid_prefix(b)
                length(ballot) <- k
                ballot
            }
        )
        theseballots <- t(theseballots)
        ballots <- rbind(ballots,theseballots)
        #min block
        theseballots <- sapply(
            seq(floor(numballots*(minority_share)*(1-minority_crossover))),
            function(x){
                if (minority_vary_majority) {
                    majorder <- sample.int(length(majority_candidates))
                }
                else {
                    majorder <- c(1:length(majority_candidates))   
                }
                if (minority_vary_minority) {
                    minorder <- sample.int(length(minority_candidates))
                }
                else {
                    minorder <- c(1:length(minority_candidates))   
                }
                chosen_type <- sample(1:nrow(ballot_type_frames$Cfirst), prob=ballot_type_frames$Cfirst$p, size=1)
                positions <- get_positions(ballot_type_frames$Cfirst,chosen_type,20)
                majpositions <- positions$majpositions[1:min(length(majority_candidates),length(positions$majpositions))]
                minpositions <- positions$minpositions[1:min(length(minority_candidates),length(positions$minpositions))]
                b <- vector()
                length(b) <- k
                b[majpositions] <- majority_candidates[majorder[1:length(majpositions)]]
                b[minpositions] <- minority_candidates[minorder[1:length(minpositions)]]
                ballot <- valid_prefix(b)
                length(ballot) <- k
                ballot
            }
        )
        theseballots <- t(theseballots)
        ballots <- rbind(ballots,theseballots)
        #min cross
        theseballots <- sapply(
            seq(floor(numballots*(minority_share)*minority_crossover)),
            function(x){
                if (minority_vary_majority) {
                    majorder <- sample.int(length(majority_candidates))
                }
                else {
                    majorder <- c(1:length(majority_candidates))   
                }
                if (minority_vary_minority) {
                    minorder <- sample.int(length(minority_candidates))
                }
                else {
                    minorder <- c(1:length(minority_candidates))   
                }
                chosen_type <- sample(1:nrow(ballot_type_frames$Wfirst), prob=ballot_type_frames$Wfirst$p, size=1)
                positions <- get_positions(ballot_type_frames$Wfirst,chosen_type,20)
                majpositions <- positions$majpositions[1:min(length(majority_candidates),length(positions$majpositions))]
                minpositions <- positions$minpositions[1:min(length(minority_candidates),length(positions$minpositions))]
                b <- vector()
                length(b) <- k
                b[majpositions] <- majority_candidates[majorder[1:length(majpositions)]]
                b[minpositions] <- minority_candidates[minorder[1:length(minpositions)]]
                ballot <- valid_prefix(b)
                length(ballot) <- k
                ballot
            }
        )
        theseballots <- t(theseballots)
        ballots <- rbind(ballots,theseballots)

        #convert to ranks e.g. 1,2,3,4
        #print(ballots)
        ranksframe = data.frame(matrix(nrow=0,ncol=k))
        names(ranksframe) <- c(minority_candidates,majority_candidates)
        for (b in c(1:nrow(ballots))){
            for (rank in which(!is.na(ballots[b,]))){
                ranksframe[b,ballots[b,rank]] = rank
            }
        }
        ranksframe
    
}

multiGroupBallotsAlternatingCrossover <- function(
    minority_crossover, majority_crossover,
    minority_candidates, majority_candidates,
    minority_vary_minority, minority_vary_majority,
    majority_vary_minority, majority_vary_majority,
    minority_share,
    numballots){
    ballots = matrix(nrow=0,ncol=length(minority_candidates)+length(majority_candidates))
    k = length(majority_candidates)+length(minority_candidates)
    colnames(ballots) <- c(1:k)
    #maj block
    theseballots <- sapply(
        seq(floor(numballots*(1-minority_share)*(1-majority_crossover))),
        function(x){
            if (majority_vary_majority) {
                majorder <- sample.int(length(majority_candidates))
            }
            else {
                majorder <- c(1:length(majority_candidates))   
            }
            if (majority_vary_minority) {
                minorder <- sample.int(length(minority_candidates))
            }
            else {
                minorder <- c(1:length(minority_candidates))   
            }
            majpositions <- c(1:length(majority_candidates))
            minpositions <- c((length(majority_candidates)+1):k)
            b <- c(1:k)
            b[majpositions] <- majority_candidates[majorder]
            b[minpositions] <- minority_candidates[minorder]
            b
        }
    )
    theseballots <- t(theseballots)
    ballots <- rbind(ballots,theseballots)
    #maj cross
    theseballots <- sapply(
        seq(floor(numballots*(1-minority_share)*majority_crossover)),
        function(x){
            if (majority_vary_majority) {
                majorder <- sample.int(length(majority_candidates))
            }
            else {
                majorder <- c(1:length(majority_candidates))   
            }
            if (majority_vary_minority) {
                minorder <- sample.int(length(minority_candidates))
            }
            else {
                minorder <- c(1:length(minority_candidates))   
            }
            if (length(minority_candidates) == length(majority_candidates)){
                majpositions <- 2*c(1:length(majority_candidates))
                minpositions <- c(2*c(1:length(majority_candidates))-1)
            }
            if (length(minority_candidates) > length(majority_candidates)){
                majpositions <- 2*c(1:length(majority_candidates))
                minpositions <- c(2*c(1:length(majority_candidates))-1, c((2*(length(majority_candidates))+1):k))
            }
            if  (length(minority_candidates) < length(majority_candidates)) {
                majpositions <- c(2*c(1:length(minority_candidates)), c((2*(length(minority_candidates))+1):k))
                minpositions <- 2*c(1:length(minority_candidates))-1
            }
            b <- c(1:k)
            b[majpositions] <- majority_candidates[majorder]
            b[minpositions] <- minority_candidates[minorder]
            b
        }
    )
    theseballots <- t(theseballots)
    ballots <- rbind(ballots,theseballots)
    #min block
    theseballots <- sapply(
        seq(floor(numballots*(minority_share)*(1-minority_crossover))),
        function(x){
            if (minority_vary_majority) {
                majorder <- sample.int(length(majority_candidates))
            }
            else {
                majorder <- c(1:length(majority_candidates))   
            }
            if (minority_vary_minority) {
                minorder <- sample.int(length(minority_candidates))
            }
            else {
                minorder <- c(1:length(minority_candidates))   
            }
            majpositions <- c((length(minority_candidates)+1):k) 
            minpositions <- c(1:length(minority_candidates))
            b <- c(1:k)
            b[majpositions] <- majority_candidates[majorder]
            b[minpositions] <- minority_candidates[minorder]
            b
        }
    )
    theseballots <- t(theseballots)
    ballots <- rbind(ballots,theseballots)
    #min cross
    theseballots <- sapply(
        seq(floor(numballots*(minority_share)*minority_crossover)),
        function(x){
            if (minority_vary_majority) {
                majorder <- sample.int(length(majority_candidates))
            }
            else {
                majorder <- c(1:length(majority_candidates))   
            }
            if (minority_vary_minority) {
                minorder <- sample.int(length(minority_candidates))
            }
            else {
                minorder <- c(1:length(minority_candidates))   
            }
            if (length(majority_candidates) == length(minority_candidates)){
                minpositions <- 2*c(1:length(minority_candidates))
                majpositions <- 2*c(1:length(minority_candidates))-1
            }
            if (length(majority_candidates) > length(minority_candidates)){
                minpositions <- 2*c(1:length(minority_candidates))
                majpositions <- c(2*c(1:length(minority_candidates))-1, c((2*(length(minority_candidates))+1):k))
            }
            if (length(majority_candidates) < length(minority_candidates)){
                minpositions <- c(2*c(1:length(majority_candidates)), c((2*(length(majority_candidates))+1):k))
                majpositions <- 2*c(1:length(majority_candidates))-1
            }
            b <- c(1:k)
            b[majpositions] <- majority_candidates[majorder]
            b[minpositions] <- minority_candidates[minorder]
            b
        }
    )
    theseballots <- t(theseballots)
    ballots <- rbind(ballots,theseballots)
    
    #convert to ranks e.g. 1,2,3,4
    #print(ballots)
    ranksframe = data.frame(matrix(nrow=0,ncol=k))
    names(ranksframe) <- c(minority_candidates,majority_candidates)
    for (b in c(1:nrow(ballots))){
        for (rank in c(1:k)){
            ranksframe[b,ballots[b,rank]] = rank
        }
    }
    ranksframe
}

#returns a list of ranked ballots using a Paired Comparison via support vector model
singleGroupBallotsPaired <- function(supportvector, candidatenames, numballots){
    ballotframe = matrix(nrow=numballots,ncol=length(candidatenames))
    #generate ballots e.g. A,B,C,D
    ballot = c(1:length(candidatenames))
    ballotsin = 1
    steps = 1
    while (ballotsin <= numballots){
        toswapleft <- sample.int(length(candidatenames)-1, size=1)
        currentcomparison <- supportvector[ballot[toswapleft]]
        proposedcomparison <- supportvector[ballot[toswapleft+1]]
        if (runif(1,0,1) < proposedcomparison/currentcomparison){
            ballot <- replace(
                ballot,
                c(toswapleft, toswapleft+1),
                ballot[c(toswapleft+1, toswapleft)]
            )
        }
        if (steps %% sample_interval && steps >= burn){
            ballotframe[ballotsin,] <- candidatenames[ballot]
            ballotsin <- ballotsin + 1
        }
        steps <-  steps+1
    }
    #convert to ranks e.g. 1,2,3,4
    ranksframe = data.frame(matrix(nrow=0,ncol=length(candidatenames)))
    names(ranksframe) <- candidatenames
    for (b in c(1:nrow(ballotframe))){
        for (rank in c(1:length(candidatenames))){
            ranksframe[b,ballotframe[b,rank]] = rank
        }
    }
    ranksframe
}

multiGroupBallotsPaired = function(supportmatrix, candidatenames, totalballots, proportions){
    ballots = data.frame(matrix(nrow=0,ncol=length(candidatenames)))
    names(ballots) <- candidatenames
    for (r in 1:length(proportions)){
        ballots = rbind(ballots,
                        singleGroupBallotsPaired(
                            supportmatrix[r,],
                            candidatenames,
                            floor(totalballots*proportions[r])
                        )
        )
    }
    ballots
}

#returns a list of ranked ballots using a Luce model
singleGroupBallots = function(supportvector, candidatenames, numballots){
    #generate ballots e.g. A,B,C,D
    ballots <- sapply(
        seq(numballots),
        function(x){sample(
            candidatenames,
            length(candidatenames),
            replace=FALSE,
            prob=supportvector
        )}
    )
    ballotframe <- t(ballots)

    #convert to ranks e.g. 1,2,3,4
    ranksframe = data.frame(matrix(nrow=0,ncol=length(candidatenames)))
    names(ranksframe) <- candidatenames
    for (b in c(1:nrow(ballotframe))){
        for (rank in c(1:length(candidatenames))){
            ranksframe[b,ballotframe[b,rank]] = rank
        }
    }
    ranksframe
}

multiGroupBallots = function(supportmatrix, candidatenames, totalballots, proportions){
    ballots = data.frame(matrix(nrow=0,ncol=length(candidatenames)))
    names(ballots) <- candidatenames
    for (r in 1:length(proportions)){
        ballots <- rbind(ballots,
                        singleGroupBallots(
                            supportmatrix[r,],
                            candidatenames,
                            floor(totalballots*proportions[r])
                        )
        )
    }
    ballots
}

rcvWinners = function(df, totalballots, proportions, seats, model_choice, ballot_type_frames=NULL, scenario_parameters=NULL){
    M <- df
    candidatenames <- names(df)
    if (model_choice == 'Luce model') {
        ranksframe <- multiGroupBallots(M, candidatenames, totalballots, proportions)
    }
    if (model_choice == 'Bradley-Terry model') {
        ranksframe <- multiGroupBallotsPaired(M, candidatenames, totalballots, proportions)
    }
    if (model_choice == 'Alternating crossover'){
        ranksframe <- multiGroupBallotsAlternatingCrossover(
            scenario_parameters$minority_crossover, scenario_parameters$majority_crossover,
            scenario_parameters$minority_candidates, scenario_parameters$majority_candidates,
            scenario_parameters$minority_vary_minority, scenario_parameters$minority_vary_majority,
            scenario_parameters$majority_vary_minority, scenario_parameters$majority_vary_majority,
            proportions[1],totalballots
        )
    }
    if (model_choice == 'Cambridge sampler'){
        ranksframe <- multiGroupBallotsCambridgeSampler(
            scenario_parameters$minority_crossover, scenario_parameters$majority_crossover,
            scenario_parameters$minority_candidates, scenario_parameters$majority_candidates,
            scenario_parameters$minority_vary_minority, scenario_parameters$minority_vary_majority,
            scenario_parameters$majority_vary_minority, scenario_parameters$majority_vary_majority,
            proportions[1],totalballots,ballot_type_frames
        )
    }
    rcvresults <- vote::stv(ranksframe, mcan = seats, eps=1, verbose=FALSE, seed=NULL)
    #write.csv(rcvresults$details, "election_rounds_6seats_6candidateseach.csv")
    pluralityelected <- plurality(ranksframe, seats)
    return(list(rcv=rcvresults$elected, plurality=pluralityelected))
}

plurality = function(ranksframe, seats){
    votematrix <- data.matrix(ranksframe)
    binaryvotematrix <- apply(votematrix, 1:2, function(x){
        if (is.na(x)) {0}
        else if (x <= seats){1}
        else {0}
    })
    binaryvotematrix = binaryvotematrix[,order(colSums(-binaryvotematrix))]
    return(colnames(binaryvotematrix)[1:seats])
}


ui <- fluidPage(theme=shinytheme("cosmo"),
    titlePanel("Ranked choice voting and representation"),

    sidebarLayout(
        sidebarPanel(
            tags$head(
                tags$style(HTML('#run{background-color:blue}'))
            ),
            numericInput("seats", "Number of seats open",3, min=0, max=10),
            numericInput("numballots","Number of ballots cast*",1000, max=10000, min=0),
            numericInput("prop1","Minority percentage",30,min=0, max=100),
            h3(""),
            "*Although 1000 ballots is low for most real-world elections, we recommend this value
            for computational efficiency and since it is high enough to simulate most effects of interest",
        ),
        
        mainPanel(
            tabsetPanel(type='tabs',
                tabPanel("Model inputs",
                         h3("Type of model"),
                         radioButtons(
                             "choice",
                             "",
                             choices=c("Luce model", "Bradley-Terry model", "Alternating crossover", "Cambridge sampler")
                         ),
                         h3("Basic inputs"),
                         uiOutput("parameterinputs"),
                         h3("Model-specific parameters"),
                         uiOutput("modelparameterinputs")
                ),
                tabPanel("Simulate a few elections",
                         h3(""),
                         actionButton("run", "Simulate an election"),
                         h3(""),
                         actionButton("reset", "Clear election history"),
                         h3("Candidates currently up for election:"),
                         verbatimTextOutput("candidatesup"),
                         h3("Winners under RCV:"),
                         verbatimTextOutput("displayrcvwinners"),
                         h3("Winners under (un-ranked) plurality:"),
                         verbatimTextOutput("displaypluralitywinners")                        
                ),
                tabPanel("Aggregate over multiple elections",
                         h3(""),
                         numericInput("numsimulations", "Number of simulated elections to run", 10, min=0, max=500),
                         actionButton("runmulti", "Run multiple elections"),
                         h3("Results:"),
                         plotOutput("displaypreferredhist"),
                         plotOutput("displayaggregates")
                ),
                tabPanel("About the models",
                         h3("Modelling voting behavior"),
                         "Voting is an incredibly complex process, and any attempt to model
                         it is bound to oversimplify the way people vote, especially
                         in ranked elections, where a voter has so many more ways to make
                         their voice heard. With this in mind, this app uses multiple models of ranking,
                         some classical and some new. They all take the same basic input.",
                         h3("Basic inputs to the models"),
                         "All the models take some basic inputs, namely the support from each group for each
                         class of candidate (minority- or majority-preferred), which should be numbers
                         between 0 and 1. These values can be estimated
                         by analyzing single-winner elections, inferred from survey results, or set to hypothetical values. 
                         Each model interprets these support values slightly differently, as detailed below.",
                         h3("How the model works: Luce model"),
                         "The Luce model is a classical model of statistical ranking
                         which has already been applied to study ranked choice voting in Ireland [1].
                         Let's suppose the support for Candidate A from a voting group is twice
                         the support for Candidate B. Under the Luce model, we take this to
                         mean that a voter is twice as likely to put Candidate A first as
                         to put Candidate B first. In fact, the model assumes that this trend
                         goes ", em("all the way down the ballot"), ". That is, if a voter
                         has not ranked either A or B yet on their ballot, then they are still
                         twice as likely to write down Candidate A's name than Candidate B when
                         they get to the next rank. See [2] for mathematical details.
                         It is unrealistic to expect the support for each class of
                         candidates to be completely uniform, so to divide up the 
                         total support among candidates of a given class, we sample from a symmetric Dirichlet distribution.
                         The Dirichlet distribution has one parameter
                         (called concentration) which can be used to choose between mostly even division of support (concentration >> 1),
                         and support concentrated on just a few candidates (concentration << 1). The default concentration value of 1.0 makes every division
                         of the support equally likely, roughly speaking. We therefore refer to this parameter as ", em("evenness of support")," in the Model input tab.",
                         h3("How the model works: Bradley-Terry"),
                         "The paired comparison model is based on the idea that what matters
                         in a ranking is who is preferred over whom. Suppose the support by a 
                         voting group for Candidate A
                         is ", em("a"), " and the support for Candidate B is ", em("b"),
                         ". This the model sets the probability that a voter in the group ranks
                         Candidate A above Candidate B at ", em("a/(a+b)"), ". 
                         The probability of a full ranking is just the product of the probabilities
                         associated with the order of each pair of candidates, with a normalizing constant
                         so that everything sums to one. See [2] for mathematical details. As with the Plackett-Luce model,
                         we use a Dirichlet sampler to divide up the support among the candidates.",
                         h3("How the model works: Alternating crossover"),
                         "This model posits that there are two kinds of voters in each group: block voters and crossover voters.
                         Block voters always vote for ingroup candidates first and then outgroup candidates. Crossover voters
                         alternate between the two classes of candidates, starting with an outgroup candidate (hence the name). This only
                         tells you the type of candidate at each position in a ballot, however, so further choices are necessary
                         to determine the exact order. For each group of voters and group of candidates, we allow two choices:
                         random ordering by each voter, or consistent ordering by each voter. This model has been used in previous analyses of
                         ranked voting by the MGGG at Tufts University (see [3,4] below).",
                         h3("How the model works: Cambridge sampler"),
                         "Each voter's first choice under this model is the same as under alternating crossover, i.e. an ingroup candidate for block voters
                         and an outgroup candidate for crossover voters. However, instead of allowing only a block ballot or alternating ballot, 
                         the types of candidates in the rest of the ballot are determined by sampling from the ballots in a decade's worth of Cambridge MA
                         ranked choice city council elections. Once the type of candidate (ingroup or outgroup) at each rank has been chosen, 
                         random ordering or consistent ordering is applied to fill in the ballot with candidate names just as with the alternating crossover model.",
                         h3("References"),
                         "[1] Gormley, I.S, and Murphy, T.B.",em("Exploring voting blocs within the Irish electorate: A mixture modeling approach."),
                         "Journal of the American Statistical Association 103.483 (2008): 1014-1027.",
                         HTML("<br>"),
                         "[2] Critchlow, D.E., Fligner, M.A. and Verducci, J.S., 1991.",
                         em("Probability models on rankings"), ". Journal of mathematical psychology, 35(3), pp.294-318.",
                         HTML("<br>"),
                         "[3] MGGG,", em("Election Reform in Lowell MA, "), a("https://mggg.org/lowell", href="https://mggg.org/lowell"),
                         HTML("<br>"),
                         "[4] MGGG,", em("Analysis of county commission elections in Yakima County WA, "), a("https://mggg.org/uploads/yakima.pdf", href="https://mggg.org/uploads/yakima.pdf")
                )
            )

        )
    )
)

server <- function(input, output, session) {
    alldisplaywinners <- reactiveValues(rcv="", plurality="", trials=0, clear=TRUE)
    candidatenames <- reactiveVal()
    electionsdone <- reactiveVal(0)
    ballot_type_frames <- reactiveValues(Wfirst = NULL, Cfirst = NULL)
    
    #force sum to one for inputs
    observeEvent(input$minmean_1,{
        updateNumericInput(session, "minmean_2", value=1-input$minmean_1)
    })
    observeEvent(input$minmean_2,{
        updateNumericInput(session, "minmean_1", value=1-input$minmean_2)
    })
    observeEvent(input$majmean_1,{
        updateNumericInput(session, "majmean_2", value=1-input$majmean_1)
    })
    observeEvent(input$majmean_2,{
        updateNumericInput(session, "majmean_1", value=1-input$majmean_2)
    })
    
    #load Cambridge ballot types
    ballot_type_frames$Wfirst <- read.csv("Wfirst_all_Cambridge_types_encoded.csv")
    ballot_type_frames$Cfirst <- read.csv("Cfirst_all_Cambridge_types_encoded.csv")
    
    #get list of candidate names
    candidatenames <- reactive({
        candidatenames <- unlist(
            lapply(seq(numcandtypes), function(x){
                sapply(seq(input[[paste('copies',x,sep='_')]]), function(y) {paste(input[[paste('name',x,sep='_')]],y,sep="")})
            }
            )
        )
    })
    
    #render input widgets for candidate types
    output$parameterinputs <- renderUI({
        alphabet <- 'ABCDEFGHIHJKLMNOPQRSTUVWXYZ'
        #all-model parameters
            tagList(
            div(
                style="display: inline-block;vertical-align:top; width: 210px;",
                textInput(paste("name", 1, sep="_"), "Minority-preferred candidate prefix", 'MIN'),
            ),
            div(style="display: inline-block;vertical-align:top; width: 40px;",HTML("<br>")),
            div(
                style="display: inline-block;vertical-align:top; width: 210px;",
                numericInput(paste("copies",1, sep="_"), "No. of minority-preferred candidates", 3, max=10),
            ),
            HTML("<br>"),
            div(
                style="display: inline-block;vertical-align:top; width: 210px;",
                numericInput(paste("minmean",1, sep="_"),"Minority support", 0.7, min=0, max=1, step=0.1),
            ),
            div(style="display: inline-block;vertical-align:top; width: 40px;",HTML("<br>")),
            div(
                style="display: inline-block;vertical-align:top; width: 210px;",
                numericInput(paste("majmean",1, sep="_"),"Majority support",0.2, min=0, max=1, step=0.1),
            ),
            hr(),
            div(
                style="display: inline-block;vertical-align:top; width: 210px;",
                textInput(paste("name", 2, sep="_"), "Majority-preferred candidate prefix", 'MAJ'),
            ),
            div(style="display: inline-block;vertical-align:top; width: 40px;",HTML("<br>")),
            div(
                style="display: inline-block;vertical-align:top; width: 210px;",
                numericInput(paste("copies",2, sep="_"), "No. of majority-preferred candidates", 3, max=10),
            ),
            HTML("<br>"),
            div(
                style="display: inline-block;vertical-align:top; width: 210px;",
                numericInput(paste("minmean",2, sep="_"),"Minority support", 0.3, min=0, max=1, step=0.1),
            ),
            div(style="display: inline-block;vertical-align:top; width: 40px;",HTML("<br>")),
            div(
                style="display: inline-block;vertical-align:top; width: 210px;",
                numericInput(paste("majmean",2, sep="_"),"Majority support",0.8, min=0, max=1, step=0.1),
            ),
            )
    })
    
    #special parameters for each model
    output$modelparameterinputs <- renderUI({
        if (input$choice == 'Luce model' || input$choice == 'Bradley-Terry model'){
            tagList(
                numericInput(paste("minstd",1, sep="_"), "Evenness of support by minority for minority-preferred candidates", 1.0, step=0.1, min=0, max=10),
                numericInput(paste("majstd",1, sep="_"), "Evenness of support by minority for majority-preferred candidates",1.0, step=0.1, min=0, max=10),
                numericInput(paste("minstd",2, sep="_"), "Evenness of support by majority for minority-preferred candidates", 1.0, step=0.1, min=0, max=10),
                numericInput(paste("majstd",2, sep="_"), "Evenness of support by majority for majority-preferred candidates",1.0, step=0.1, min=0, max=10),
            )
        }
        else {
            tagList(
                checkboxInput("minority_order_minority", "Minority voters agree on the order of the minority candidates.", TRUE),
                checkboxInput("minority_order_majority", "Minority voters agree on the order of the majority candidates.", TRUE),
                checkboxInput("majority_order_minority", "Majority voters agree on the order of the minority candidates.", TRUE),
                checkboxInput("majority_order_majority", "Majority voters agree on the order of the majority candidates.", TRUE)
                
            )
        }
    })

    #gets rcvwinners for a single election when button pressed
    winners <- eventReactive(input$run, {
        alldisplaywinners$trials <- alldisplaywinners$trials + 1
        #build dataframe for parameters
        cands <- candidatenames()
        df <- data.frame(matrix(ncol=length(cands),nrow=2))
        names(df) <- cands
        
        #------models that use support vectors------
        if (input$choice == 'Luce model' || input$choice == 'Bradley-Terry model') {
            #sample support values
            minsupportvaluesperturb <- unlist(
                sapply(seq(numcandtypes), function(x){
                    rdirichlet(1, rep(input[[paste('minstd',x,sep='_')]], input[[paste('copies',x,sep='_')]]))})
            )
            minsupportvaluesmean <- unlist(
                sapply(seq(numcandtypes), function(x){
                    rep(input[[paste('minmean',x,sep='_')]]/input[[paste('copies',x,sep='_')]], input[[paste('copies',x,sep='_')]])})
            )
            minsupportvalues <- minsupportvaluesmean*minsupportvaluesperturb
            majsupportvaluesperturb <- unlist(
                sapply(seq(numcandtypes), function(x){
                    rdirichlet(1, rep(input[[paste('minstd',x,sep='_')]], input[[paste('copies',x,sep='_')]]))})
            )
            majsupportvaluesmean <- unlist(
                sapply(seq(numcandtypes), function(x){
                    rep(input[[paste('majmean',x,sep='_')]]/input[[paste('copies',x,sep='_')]], input[[paste('copies',x,sep='_')]])})
            )
            majsupportvalues <- majsupportvaluesmean*majsupportvaluesperturb
            df[1,] = minsupportvalues
            df[2,] = majsupportvalues
            
            #compute RCV winners
            #small corrections to zero values
            for (r in 1:nrow(df)){
                for (c in 1:ncol(df)){
                    if (df[r,c] == 0){
                        df[r,c] <- 1e-7
                    }
                }
            }
            #normalize with row sums
            for (r in 1:nrow(df)){
                rowsum <- sum(df[r,])
                for (c in 1:ncol(df)){
                    df[r,c] <- df[r,c]/rowsum
                }
            }
            #check for more seats than candidates
            if (input$seats > ncol(df)){
                allwinners <- list(
                    rcv="Error: more seats than candidates.",
                    plurality="Error: more seats than candidates."
                )
                return(allwinners)
            }
            #compute winners
            allwinners <- rcvWinners(
                df,
                input$numballots,
                c(input$prop1/100,1-input$prop1/100),
                input$seats,
                input$choice
            )
        }
        
        #----------models that use scenarios------------------------
        if (input$choice == 'Alternating crossover' || input$choice == 'Cambridge sampler'){
            scenario_parameters <- list(
                minority_crossover=input$minmean_2, majority_crossover=input$majmean_1,
                minority_candidates=candidatenames()[1:input$copies_1], majority_candidates=candidatenames()[(input$copies_1+1):length(candidatenames())],
                minority_vary_minority=!(input$minority_order_minority), minority_vary_majority=!(input$minority_order_majority),
                majority_vary_minority=!(input$majority_order_minority), majority_vary_majority=!(input$majority_order_majority)
            ) 
            allwinners <- rcvWinners(
                NULL,
                input$numballots,
                c(input$prop1/100,1-input$prop1/100),
                input$seats,
                input$choice,
                scenario_parameters=scenario_parameters,
                ballot_type_frames=ballot_type_frames
            )
        }
        
        return(allwinners)
    })
    
    #run multiple elections and get aggregate data
    aggregateplots <- eventReactive(input$runmulti, {
        #build dataframe for parameters
        cands <- candidatenames()
        df <- data.frame(matrix(ncol=length(cands),nrow=2))
        names(df) <- cands
        candtypevalues <- unlist(
            lapply(seq(numcandtypes), function(x){
                rep(globalcandtypes[[paste('preferred',x,sep='_')]],input[[paste('copies',x,sep='_')]])})
        )
        #check for more seats than candidates
        if (input$seats > ncol(df)){
            return(plot(main='Error: more seats than candidates'))
        }
        
        #compute winners
        winnerdf <- data.frame(matrix(0,ncol=length(cands),nrow=1))
        minorityseats <- c()
        names(winnerdf) <- c(cands)
        electionsdone(0)
        withProgress(message='Simulating elections...',value=0,{
        for (sim in seq(input$numsimulations)){
            #-------models that use support vectors----
            if (input$choice == 'Luce model' || input$choice == 'Bradley-Terry model'){
                #sample support values
                minsupportvaluesperturb <- unlist(
                    sapply(seq(numcandtypes), function(x){
                        rdirichlet(1, rep(input[[paste('minstd',x,sep='_')]], input[[paste('copies',x,sep='_')]]))})
                )
                minsupportvaluesmean <- unlist(
                    sapply(seq(numcandtypes), function(x){
                        rep(input[[paste('minmean',x,sep='_')]]/input[[paste('copies',x,sep='_')]], input[[paste('copies',x,sep='_')]])})
                )
                minsupportvalues <- minsupportvaluesmean*minsupportvaluesperturb
                majsupportvaluesperturb <- unlist(
                    sapply(seq(numcandtypes), function(x){
                        rdirichlet(1, rep(input[[paste('minstd',x,sep='_')]], input[[paste('copies',x,sep='_')]]))})
                )
                majsupportvaluesmean <- unlist(
                    sapply(seq(numcandtypes), function(x){
                        rep(input[[paste('majmean',x,sep='_')]]/input[[paste('copies',x,sep='_')]], input[[paste('copies',x,sep='_')]])})
                )
                majsupportvalues <- majsupportvaluesmean*majsupportvaluesperturb
                df[1,] = minsupportvalues
                df[2,] = majsupportvalues
                
                #small corrections to zero values
                for (r in 1:nrow(df)){
                    for (c in 1:ncol(df)){
                        if (df[r,c] == 0){
                            df[r,c] <- 1e-7
                        }
                    }
                }
                
                #normalize with row sums
                for (r in 1:nrow(df)){
                    rowsum <- sum(df[r,])
                    for (c in 1:ncol(df)){
                        df[r,c] <- df[r,c]/rowsum
                    }
                }
                
                #compute winners and store
                allwinners <- rcvWinners(
                    df,
                    input$numballots,
                    c(input$prop1/100,1-input$prop1/100),
                    input$seats,
                    input$choice
                )
                rcvwinners <- allwinners$rcv
            }
            #----------models that use scenarios------------------------
            if (input$choice == 'Alternating crossover' || input$choice == 'Cambridge sampler'){
                scenario_parameters <- list(
                    minority_crossover=input$minmean_2, majority_crossover=input$majmean_1,
                    minority_candidates=candidatenames()[1:input$copies_1], majority_candidates=candidatenames()[(input$copies_1+1):length(candidatenames())],
                    minority_vary_minority=!(input$minority_order_minority), minority_vary_majority=!(input$minority_order_majority),
                    majority_vary_minority=!(input$majority_order_minority), majority_vary_majority=!(input$majority_order_majority)
                ) 
                allwinners <- rcvWinners(
                    NULL,
                    input$numballots,
                    c(input$prop1/100,1-input$prop1/100),
                    input$seats,
                    input$choice,
                    scenario_parameters=scenario_parameters,
                    ballot_type_frames=ballot_type_frames
                )
                rcvwinners <- allwinners$rcv
            }
            
            #-----count minority winners-----
            minoritywinners <- 0
            for (i in seq(length(cands))){
                if (cands[i] %in% rcvwinners){
                    winnerdf[1,cands[i]] <- winnerdf[1,cands[i]]+1
                    if (candtypevalues[i]){
                        minoritywinners <- minoritywinners + 1
                    }
                }
            }
            minorityseats <- c(minorityseats, minoritywinners)
            incProgress(1/input$numsimulations, detail=paste("",sim))
        }
        }
        )
        return(list(bar=winnerdf,hist=minorityseats))
    })
    
    output$displayaggregates <- renderPlot({
        data <- aggregateplots()
        barplot(
                sapply(seq(ncol(data$bar)),function(x){(data$bar)[1,x]}),
                names.arg=names(data$bar),
                main='Election results by candidate',
                ylab='# times elected'
        )
    })
    
    output$displaypreferredhist <- renderPlot({
        data <- aggregateplots()
        hist(data$hist, col='lightblue',xlab='# minority-preferred candidates elected',
             main='Minority-preferred candidates elected', ylab='# elections',
             breaks=seq(0,input$seats+1)-0.5)
    })
    
    #rcv display code
    output$displayrcvwinners <- reactive({
        alldisplaywinners$rcv
    })
    observeEvent(input$run, {
        elected <- winners()
        displaywinners = alldisplaywinners$rcv
        displaywinners <- paste(displaywinners, "Election #", alldisplaywinners$trials, ": ", sep="")
        for (w in 1:(length(elected$rcv)-1)){
            displaywinners = paste(displaywinners, elected$rcv[w], ", ", sep="")
        }
        displaywinners = paste(displaywinners, elected$rcv[length(elected$rcv)], "\n", sep="")
        alldisplaywinners$rcv <- displaywinners
    })

    #plurality display code
    output$displaypluralitywinners <- reactive({
        alldisplaywinners$plurality
    })
    observeEvent(input$run, {
        elected <- winners()
        displaywinners = alldisplaywinners$plurality
        displaywinners <- paste(displaywinners, "Election #", alldisplaywinners$trials, ": ", sep="")
        for (w in 1:(length(elected$plurality)-1)){
            displaywinners <- paste(displaywinners, elected$plurality[w], ", ", sep="")
        }
        displaywinners = paste(displaywinners, elected$plurality[length(elected$plurality)], "\n", sep="")
        alldisplaywinners$plurality <- displaywinners
    })
    
    #reset
    observeEvent(input$reset, {
        alldisplaywinners$rcv <- ""
        alldisplaywinners$plurality <- ""
        alldisplaywinners$trials <- 0
    })
    
    output$candidatesup <- eventReactive(input$run, {
        displaycandidates = ""
        for (n in 1:(length(candidatenames())-1)){
            displaycandidates = paste(displaycandidates,candidatenames()[n],", ", sep="")
        }
        displaycandidates = paste(displaycandidates,candidatenames()[length(candidatenames())])
        displaycandidates
    })
}

shinyApp(ui = ui, server = server)
