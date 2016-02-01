source("read.admb.R")
run.Simulation=function(N=10)
{
    theta <<- NULL
    for(i in sample(1:1000,N))
    {
        arg = paste("./pm -sim", i)
        # arg = paste("./pm -sim", i)
        system(arg)
        print(arg)
        P=read.fit("pm")
        theta<<-rbind(theta, P$est[1:5])
    }

}
run.Simulation()
tvalues=c(0.35,2000,0.001,sqrt(1/25),sqrt(1/25))
names=c("r","K","q",expression(sigma),expression(tau))
trtheta=cbind(theta[,1],exp(theta[,2:3]),sqrt(1/exp(theta[,4:5])))
bias=trtheta*0
for(i in 1:5)
{
    bias[,i]=(trtheta[,i]-tvalues[i])/trtheta[,i]*100
}
    
boxplot(bias,names=names)

