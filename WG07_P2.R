# Group 7 : (Jiahe Sun) (Lingda Wang) (Qi Yu)
# https://github.com/Linda0628/SP_practical2.git

COVID19_MODEL <- function(  
          # initial susceptible population size
          N=5.5*1e6, 
          # the total day
          TIME_SPAN=100, 
          # the number of initial exposed people
          INIT_EXPOSED_NUM=10, 
          # overall viral infectivity parameter
          LAMBDA=0.4/N, 
          #10% of the population with the lowest βi values
          CAUTIOUS_RATE=0.1, 
          #random sample of 0.1% of the population
          SAMPLE_RARE=1e-3, 
          # Once exposed they have a daily probability of 1/3 of entering the infectious state
          E2I = 1/3, 
          # The infectious have a daily probability of 1/5 of leaving the infectious state
          # # (recovery or transition to serious disease)
          I2R = 1/5) 

{
  # generate beta_i
  BETA <- rlnorm(N, 0, 0.5)
  BETA <- BETA/mean(BETA)
  
  ########## Records ##########
  # new infections each day
  new_nums <- rep(0, TIME_SPAN)
  # the number of new infections among cautious people with the lowest beta_i values
  new_nums_cautious  <- rep(0, TIME_SPAN)
  # the number of new infections among sampled SAMPLE_RARE people
  new_nums_sampled  <- rep(0, TIME_SPAN)
  
  # mark the top-10% cautious people
  cautious_ids <- order(BETA,decreasing=FALSE)[1:(CAUTIOUS_RATE*N)]
  
  # mark 0.1% sample
  sample_ids <- sample(1:N, N*SAMPLE_RARE, replace = F)
  
  # Start the epidemic by setting 10 randomly chosen people to the exposed (E) state
  E_ids <- sample(1:N, INIT_EXPOSED_NUM)
  
  ## initialize to susceptible state
  x <- rep(0,N) 
  ## create some E state
  x[E_ids] <- 1 
  ## set up storage for pop in each state 
  S <- E <- I <- R <- rep(0,TIME_SPAN) 
  ## initialize
  S[1] <- N-INIT_EXPOSED_NUM ; E[1] <- INIT_EXPOSED_NUM  
  
  for (i in 2:TIME_SPAN) {
    sum_beta = sum(BETA[x==2])
    infected_prob  <- rep(0,N)
    infected_prob[which(x==0)] = LAMBDA * sum_beta * BETA[x==0]
    
    u = runif(N)
    # the index in the S stage (before) one day Infection process
    index1 <- which(x==0) 
    ## I -> R with prob delta
    x[x==2&u<I2R] <- 3 
    ## E -> I with prob gamma 
    x[x==1&u<E2I] <- 2 
    ## S -> E with prob beta*I[i-1]
    x[x==0&u<infected_prob] <- 1
    # the index in the S stage (after) one day Infection process
    index2 <- which(x==0) 
    # the index that from S stage to E stage (after) one day Infection process
    index <- index1[index1 %in% index2 == F] 
    
    # the number of people in each stage (after) one day Infection process
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
    
    #(1) new whole population
    new_nums[i] <- S[i-1] - S[i]
    
    #(2) new cautious data (10%)
    index_cautious = cautious_ids[cautious_ids %in% index]
    new_nums_cautious[i] = length(index_cautious)
    
    #(3) new random sample data (0.1%)
    index_sample = sample_ids[sample_ids %in% index]
    new_nums_sampled[i] = length(index_sample)
  }
  
  list(new_nums=new_nums,
       new_nums_cautious=new_nums_cautious,
       new_nums_sampled=new_nums_sampled)
}


par(mfrow=c(2,2),mar=c(2.7,2.7,2,1))
i=1 ; runs = 10
while(i<=runs){
  res <- COVID19_MODEL() 
  
  # standardized three different data
  new_nums = scale(res$new_nums)
  new_nums_cautious = scale(res$new_nums_cautious)
  new_nums_sampled = scale(res$new_nums_sampled)
  
  # The maximum of three different data
  y_max=max(max(new_nums),max(new_nums_cautious),max(new_nums_sampled))
  # plot new whole population distribution
  plot(new_nums,ylim=c(-1,y_max),xlim=c(0,110),
       type = 'l', lwd=2, col=1,
       main = c('Standardized New Infections Daily')) 
  # Add new cautious data (10%) distribution
  lines(new_nums_cautious,col=2,lty=1,lwd=2) 
  # Add new random sample data (0.1%) distribution
  lines(new_nums_sampled,col=3,lty=1,lwd=2) 
  
  # peaks for three different curves
  abline(v=c(which.max(new_nums),
             which.max(new_nums_cautious), 
             which.max(new_nums_sampled)),
         col=c(1,2,3),lty=c(2,2,2))
  abline(h= c(max(new_nums),
              max(new_nums_cautious),
              max(new_nums_sampled)),
         col=c(1,2,3),lty=c(2,2,2))
  
  # mark the peaks' coordinate in the plot
  text(c(which.max(new_nums),
         which.max(new_nums_cautious),
         which.max(new_nums_sampled)),
       c(max(new_nums),
         max(new_nums_cautious),
         max(new_nums_sampled)),
       c(paste("(" , which.max(new_nums) , "," , max(new_nums) , ")"),
         paste("(" , which.max(new_nums_cautious) , "," , max(new_nums_cautious) , ")"),
         paste("(" , which.max(new_nums_sampled) , "," , max(new_nums_sampled) , ")")),
       col=c(1,2,3),cex=c(0.7,0.7,0.7))
  # label
  legend('center',c("whole population","cautious 10%","0.1% random sample"),
         col=c(1,2,3),bty='n',lty=c(1,1,1))
  i=i+1
}

################ In the picture:
## 1. The red curve(using zoe app) grows more slowly than the black curve( in the first 90 days). 
##   However, this group over a longer period of time 
##   and at a higher peak than the general population. 

## 2. It is obvious that the green curve(random sample) is basically consistent with the black curve(population), 
##    indicating that the normalized random sample is correct. 

## 3. From the 60th day to the 90th day, there was a significant increase in the curve

## 4. Compare the pictures given in the question, The figures we created in SEIR model coincide closely
##    with the result of the actual daily Covid deaths data in the left graph
##    (Among all 10 graphs, the first graph in 'visualization result 1.png' is more fitted)

## 5. Combined with the actual situation diagram, we can see that in the 0-100 day , the peak time is about 80-90 day.
##    this can be explained by R numbers(it was greater than 1 during this period),which means the 
##    infected people have high probability to transmit virus.
##    After that, the government imposed a national lockdown on day 95, R number drops below 1.
##    As a result, the number of new daily cases has dropped sharply.



################ Inference:
## 1. As can be seen from the figure, people with a strong sense of protection may not be able to 
##    keep away from the virus, but may take longer to spread the disease for some reason 
##    (they miss the time of herd immunity, making their immunity lower than the population).

## 2. We were surprised to find that people with stronger awareness of protection had a higher peak value. 
##    One possibility is that the model selected cannot reflect the actual situation, 
##    or the transfer probability may not be constant and may increase or decrease with the development of the epidemic.
##    (In fact ,from the change of R number in actual scenarios, it indicates that E2I is not a constant) 
##    While log Normal distribution may not fit the description of contact rate, we may try another model.



